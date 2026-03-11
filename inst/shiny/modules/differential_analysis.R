# =====================================================
# 差异分析模块
# =====================================================

differential_analysis_server <- function(input, output, session) {

  # 获取数据输入模块的函数
  data_input <- data_input_server(input, output, session)

  # --- 样本数量显示 ---
  output$sample_count_display <- renderUI({
    req(data_input$raw_data(), input$control_group, input$treat_group)

    # 获取组别信息
    groups_list <- list(
      Control = input$control_group,
      Treatment = input$treat_group
    )

    # 计算每组样本数
    group_counts <- sapply(groups_list, length)
    min_replicates <- min(group_counts)

    # 确定分析方法
    if (min_replicates >= 3) {
      method_text <- "limma-voom (样本充足)"
      method_color <- "success"
    } else {
      method_text <- "edgeR (样本较少)"
      method_color <- "warning"
    }

    tagList(
      tags$div(class = "alert alert-info",
               tags$h5("📊 样本统计与对比设置"),
               tags$p(tags$strong("对照组 (Control): "), group_counts["Control"], " 个样本"),
               tags$p(tags$strong("处理组 (Treatment): "), group_counts["Treatment"], " 个样本"),
               tags$p(tags$strong("最小重复数: "), min_replicates),
               tags$hr(),
               tags$p(tags$strong("分析方法: "),
                      tags$span(class = paste0("text-", method_color), method_text)),
               tags$p(class = "text-muted small",
                      "规则: 每组样本数≥3时使用limma-voom，<3时使用edgeR"),
               tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
               tags$p(tags$strong("🔄 对比方向: "),
                      tags$span(class = "text-primary", "Treatment vs Control"),
                      tags$br(),
                      tags$small(class = "text-muted",
                                "log2FC > 0: 基因在处理组中上调表达",
                                tags$br(),
                                "log2FC < 0: 基因在处理组中下调表达")
               )
      )
    )
  })

  # --- 差异分析函数 ---
  perform_differential_analysis <- function(df_use, group, min_replicates) {
    # 数据预处理
    dge <- DGEList(counts = df_use, group = group)
    dge <- calcNormFactors(dge)

    # 过滤低表达基因
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]

    if (min_replicates >= 3) {
      # limma-voom 分析流程
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)

      v <- voom(dge, design, plot = FALSE)
      fit <- lmFit(v, design)

      # 设置对比（使用字符串引用design矩阵的列名）
      cm <- makeContrasts(
        Treatment_vs_Control = Treatment - Control,
        levels = design
      )

      fit2 <- contrasts.fit(fit, cm)
      fit2 <- eBayes(fit2)

      # 获取结果
      res <- topTable(fit2, coef = "Treatment_vs_Control", number = Inf)
      res$GeneID <- rownames(res)
      res <- res %>%
        dplyr::rename(
          log2FoldChange = logFC,
          pvalue = P.Value,       # limma的原始p值
          padj = adj.P.Val,       # limma的BH校正p值
          t_stat = t              # moderated t统计量
        )

      # 验证关键列是否存在
      required_cols <- c("log2FoldChange", "pvalue", "padj")
      missing_cols <- setdiff(required_cols, colnames(res))
      if (length(missing_cols) > 0) {
        stop(sprintf("limma-voom结果缺少必要列: %s。现有列: %s",
                     paste(missing_cols, collapse = ", "),
                     paste(colnames(res), collapse = ", ")))
      }

    } else {
      # edgeR 分析流程
      # 明确指定对比方向：Treatment vs Control
      # edgeR的exactTest默认对比是第2个水平 vs 第1个水平
      # 由于我们设置了levels = c("Control", "Treatment")，所以默认就是Treatment vs Control
      if (min_replicates > 1) {
        # 有重复时估计离散度
        dge <- estimateDisp(dge)
        # 明确指定对比：Treatment vs Control
        et <- exactTest(dge, pair = c("Control", "Treatment"))
      } else {
        # 无重复时使用固定离散度
        user_disp_sqrt <- 0.1  # 可以通过UI设置
        # 明确指定对比：Treatment vs Control
        et <- exactTest(dge, pair = c("Control", "Treatment"), dispersion = user_disp_sqrt^2)
      }

      # 获取结果
      res <- topTags(et, n = Inf)$table
      res$GeneID <- rownames(res)
      res <- res %>%
        dplyr::rename(
          log2FoldChange = logFC,
          pvalue = PValue,
          padj = FDR  # edgeR的topTags返回FDR列（Benjamini-Hochberg校正）
        )

      # 验证关键列是否存在
      required_cols <- c("log2FoldChange", "pvalue", "padj")
      missing_cols <- setdiff(required_cols, colnames(res))
      if (length(missing_cols) > 0) {
        stop(sprintf("edgeR结果缺少必要列: %s。现有列: %s",
                     paste(missing_cols, collapse = ", "),
                     paste(colnames(res), collapse = ", ")))
      }
    }

    # 添加基础均值（使用标准化后的CPM值）
    res$baseMean <- rowMeans(edgeR::cpm(dge, log = FALSE, prior.count = 1))

    # 添加logCPM
    res$logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1) %>%
      rowMeans()

    return(res)
  }

  # --- 执行差异分析 ---
  deg_results <- eventReactive(input$analyze, {
    req(data_input$raw_data(), input$control_group, input$treat_group)

    # 获取数据
    df <- data_input$raw_data()
    ctrl <- input$control_group
    trt <- input$treat_group

    # 验证输入
    if (length(ctrl) == 0 || length(trt) == 0) {
      showNotification("请至少选择一个对照组和处理组样本", type = "error")
      return(NULL)
    }

    # 检查样本重叠
    if (length(intersect(ctrl, trt)) > 0) {
      showNotification("对照组和处理组不能有重叠样本", type = "error")
      return(NULL)
    }

    # 准备数据
    df_use <- df[, c(ctrl, trt)]
    group <- factor(c(rep("Control", length(ctrl)),
                  rep("Treatment", length(trt))),
                levels = c("Control", "Treatment"))

    min_replicates <- min(length(ctrl), length(trt))

    # 执行差异分析
    tryCatch({
      res <- perform_differential_analysis(df_use, group, min_replicates)

      # === LogFC方向验证和提示 ===
      # 计算上调和下调基因数量，用于方向验证
      n_up <- sum(res$log2FoldChange > 0, na.rm = TRUE)
      n_down <- sum(res$log2FoldChange < 0, na.rm = TRUE)
      n_total <- n_up + n_down

      # 计算显著性基因数量
      n_significant <- sum(res$padj < input$pval_cutoff & abs(res$log2FoldChange) > input$log2fc_cutoff, na.rm = TRUE)
      n_up_sig <- sum(res$padj < input$pval_cutoff & res$log2FoldChange > input$log2fc_cutoff, na.rm = TRUE)
      n_down_sig <- sum(res$padj < input$pval_cutoff & res$log2FoldChange < -input$log2fc_cutoff, na.rm = TRUE)

      # p值统计
      pval_stats <- summary(res$pvalue)
      padj_stats <- summary(res$padj)

      # 显示对比方向和统计信息
      cat("\n========== 差异分析结果摘要 ==========\n")
      cat("分析方法:", if(min_replicates >= 3) "limma-voom" else "edgeR", "\n")
      cat("对比组别: Treatment vs Control\n")
      cat("含义: 相对于Control组，Treatment组的基因表达变化\n")
      cat("----------------------------------------\n")
      cat("P值校正方法: Benjamini-Hochberg (BH FDR)\n")
      cat("筛选阈值: padj <", input$pval_cutoff, "且 |log2FC| >", input$log2fc_cutoff, "\n")
      cat("----------------------------------------\n")
      cat("总体分布:\n")
      cat(sprintf("  上调基因 (log2FC > 0): %d (%.1f%%)\n", n_up, 100*n_up/n_total))
      cat(sprintf("  下调基因 (log2FC < 0): %d (%.1f%%)\n", n_down, 100*n_down/n_total))
      cat("----------------------------------------\n")
      cat("显著差异基因:\n")
      cat(sprintf("  总计: %d\n", n_significant))
      cat(sprintf("  上调: %d\n", n_up_sig))
      cat(sprintf("  下调: %d\n", n_down_sig))
      cat("----------------------------------------\n")
      cat("P值分布:\n")
      cat(sprintf("  最小值: %.2e\n", pval_stats["Min."]))
      cat(sprintf("  中位数: %.2e\n", pval_stats["Median"]))
      cat(sprintf("  最大值: %.2e\n", pval_stats["Max."]))
      cat("校正P值 (FDR) 分布:\n")
      cat(sprintf("  最小值: %.2e\n", padj_stats["Min."]))
      cat(sprintf("  中位数: %.2e\n", padj_stats["Median"]))
      cat(sprintf("  最大值: %.2e\n", padj_stats["Max."]))
      cat("========================================\n\n")

      # 添加差异状态
      res$Status <- ifelse(
        res$padj < input$pval_cutoff & abs(res$log2FoldChange) > input$log2fc_cutoff,
        ifelse(res$log2FoldChange > 0, "Up", "Down"),
        "Not DE"
      )

      # 添加t统计量（如果没有从topTable获取到）
      if (!"t_stat" %in% colnames(res)) {
        # 对于edgeR结果，使用近似t统计量
        res$t_stat <- qnorm(1 - res$pvalue/2) * sign(res$log2FoldChange)
      }

      # 基因注释
      anno <- data_input$annotate_genes(res$GeneID, input$species_select)

      if (!is.null(anno)) {
        # 改进的GeneID清理逻辑 - 保留Ensembl ID的版本号
        clean_geneid <- res$GeneID
        clean_geneid <- trimws(clean_geneid)
        clean_geneid <- gsub("[\t\n\r]", "", clean_geneid)
        # 对于非Ensembl ID，才移除特殊字符
        non_ensembl <- !grepl("^ENS", clean_geneid, ignore.case = TRUE)
        clean_geneid[non_ensembl] <- gsub("[^[:alnum:]]", "", clean_geneid[non_ensembl])

        # 清理anno中的列名
        anno_clean <- anno
        if ("SYMBOL" %in% colnames(anno_clean)) {
          # 使用与原始数据相同的清理逻辑
          anno_clean$SYMBOL_CLEAN <- gsub("[^[:alnum:]]", "", anno_clean$SYMBOL)
          if (input$species_select == "Mm") {
            anno_clean$SYMBOL_CLEAN <- sapply(anno_clean$SYMBOL_CLEAN, function(x) {
              if (grepl("^[A-Za-z]", x) && nchar(x) > 0) {
                paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
              } else {
                x
              }
            }, USE.NAMES = FALSE)
          } else {
            anno_clean$SYMBOL_CLEAN <- toupper(anno_clean$SYMBOL_CLEAN)
          }
        }

        # 初始化结果列
        if (!"SYMBOL" %in% colnames(res)) res$SYMBOL <- res$GeneID
        if (!"ENTREZID" %in% colnames(res)) res$ENTREZID <- NA

        # 🔥 修复：记录原始GeneID，用于后续显示
        res$Original_GeneID <- res$GeneID

        # 🔥 新的匹配策略：优先使用ENSEMBL列匹配Ensembl ID
        cat("开始SYMBOL匹配流程...\n")
        cat("注释数据列名:", paste(colnames(anno_clean), collapse=", "), "\n")

        if ("ENSEMBL" %in% colnames(anno_clean)) {
          # 检查哪些基因是Ensembl ID
          is_ensembl_id <- grepl("^ENS", clean_geneid, ignore.case = TRUE)
          n_ensembl <- sum(is_ensembl_id)
          cat("发现", n_ensembl, "个Ensembl ID格式的基因\n")

          if (any(is_ensembl_id)) {
            # 调试：显示前5个Ensembl ID
            ensembl_ids <- clean_geneid[is_ensembl_id]
            cat("前5个Ensembl ID:", paste(head(ensembl_ids, 5), collapse=", "), "\n")
            cat("注释数据库中ENSEMBL列数量:", length(anno_clean$ENSEMBL), "\n")

            # 直接通过ENSEMBL列匹配
            ensembl_match_idx <- match(clean_geneid[is_ensembl_id], anno_clean$ENSEMBL)
            matched_ensembl <- !is.na(ensembl_match_idx)
            n_matched <- sum(matched_ensembl)

            cat("通过ENSEMBL列匹配成功:", n_matched, "/", n_ensembl, "个基因\n")

            if (any(matched_ensembl)) {
              ensembl_indices <- which(is_ensembl_id)
              indices_to_update <- ensembl_indices[matched_ensembl]
              res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[ensembl_match_idx[matched_ensembl]]
              res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[ensembl_match_idx[matched_ensembl]]
            }

            # 🔥 调试：显示未匹配的基因
            if (n_matched < n_ensembl) {
              unmatched_ids <- clean_geneid[is_ensembl_id][!matched_ensembl]
              cat("警告:", sum(!matched_ensembl), "个Ensembl ID未能匹配\n")
              cat("未匹配示例（前5个）:", paste(head(unmatched_ids, 5), collapse=", "), "\n")
            }
          }
        } else {
          cat("错误：注释数据中没有ENSEMBL列！\n")
        }


        # 第一步：尝试SYMBOL匹配（对于非Ensembl ID或未匹配的基因）
        if ("SYMBOL" %in% colnames(anno_clean)) {
          # 只对尚未匹配ENTREZID的基因尝试SYMBOL匹配
          unmatched <- is.na(res$ENTREZID) | res$ENTREZID == ""

          if (any(unmatched)) {
            match_idx <- match(clean_geneid[unmatched], anno_clean$SYMBOL_CLEAN)
            matched_genes <- !is.na(match_idx)

            if (any(matched_genes)) {
              unmatched_indices <- which(unmatched)
              indices_to_update <- unmatched_indices[matched_genes]
              res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[match_idx[matched_genes]]
              res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[match_idx[matched_genes]]
              cat("通过SYMBOL列匹配成功:", sum(matched_genes), "个基因\n")
            }
          }

          # 🔥 第二步：对于仍未匹配的基因，尝试使用ENTREZID反向查询
          still_unmatched <- is.na(res$ENTREZID) | res$ENTREZID == ""
          if (any(still_unmatched)) {
            # 尝试通过ENTREZID匹配
            entrez_match_idx <- match(res$GeneID[still_unmatched], anno_clean$ENTREZID)
            matched_entrez <- !is.na(entrez_match_idx)

            if (any(matched_entrez)) {
              unmatched_indices <- which(still_unmatched)
              indices_to_update <- unmatched_indices[matched_entrez]
              res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[entrez_match_idx[matched_entrez]]
              res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[entrez_match_idx[matched_entrez]]
              cat("通过ENTREZID反向匹配成功:", sum(matched_entrez), "个基因\n")
            }
          }
        }

      }

      # 确保有SYMBOL列
      if (!"SYMBOL" %in% colnames(res)) res$SYMBOL <- res$GeneID
      if (!"ENTREZID" %in% colnames(res)) res$ENTREZID <- NA

      # 过滤假基因
      res <- data_input$filter_pseudo_genes(res)

      # 改进的基因去重逻辑 - 保留统计显著性最高的基因
      if (any(duplicated(res$SYMBOL))) {
        # 按照p值和log2FoldChange的显著性排序
        res <- res %>%
          dplyr::arrange(SYMBOL, padj, abs(log2FoldChange)) %>%
          dplyr::distinct(SYMBOL, .keep_all = TRUE)

        # 记录去重信息
        n_duplicates <- sum(duplicated(res$SYMBOL))
        if (n_duplicates > 0) {
          cat(sprintf("移除了 %d 个重复的基因记录\n", n_duplicates))
        }
      }

      return(res)

    }, error = function(e) {
      showNotification(paste("差异分析失败:", e$message), type = "error")
      return(NULL)
    })
  })

  # 增强的列映射函数
  enhanced_column_mapping <- function(df) {
    cat("检查上传的差异基因文件列结构...\n")
    cat("原始列名:", paste(colnames(df), collapse = ", "), "\n")

    # 可能的列名映射
    column_mappings <- list(
      log2FoldChange = c("log2FoldChange", "log2FC", "avg_log2FC", "logFC", "log2_fold_change", "log2fc", "log2fc_adj"),
      pvalue = c("pvalue", "p_val", "p.value", "P.Value", "pvalue_adj"),
      padj = c("padj", "p_val_adj", "p_adj", "adj.P.Val", "pvalue_adj", "FDR"),
      GeneID = c("GeneID", "gene", "Gene", "SYMBOL", "symbol", "gene_symbol", "ensembl", "ENSEMBL")
    )

    # 检查并重命名列
    for (target_col in names(column_mappings)) {
      possible_names <- column_mappings[[target_col]]
      found <- FALSE

      for (col_name in possible_names) {
        if (col_name %in% colnames(df)) {
          if (col_name != target_col) {
            cat("  重命名列:", col_name, "->", target_col, "\n")
            colnames(df)[colnames(df) == col_name] <- target_col
          } else {
            cat("  找到列:", target_col, "\n")
          }
          found <- TRUE
          break
        }
      }

      if (!found) {
        cat("  ⚠️  缺失列:", target_col, "\n")
      }
    }

    # 确保log2FoldChange是数值类型
    if ("log2FoldChange" %in% colnames(df)) {
      if (!is.numeric(df$log2FoldChange)) {
        cat("  转换log2FoldChange为数值类型\n")
        original_type <- class(df$log2FoldChange)[1]
        df$log2FoldChange <- as.numeric(as.character(df$log2FoldChange))
        n_na <- sum(is.na(df$log2FoldChange))
        if (n_na > 0) {
          warning(sprintf("log2FoldChange列从%s转换为数值时产生了%d个NA值", original_type, n_na))
        }
      }
    }

    # 确保pvalue和padj是数值类型
    for (col in c("pvalue", "padj")) {
      if (col %in% colnames(df)) {
        if (!is.numeric(df[[col]])) {
          cat("  转换", col, "为数值类型\n")
          original_type <- class(df[[col]])[1]
          df[[col]] <- as.numeric(as.character(df[[col]]))
          n_na <- sum(is.na(df[[col]]))
          if (n_na > 0) {
            warning(sprintf("%s列从%s转换为数值时产生了%d个NA值", col, original_type, n_na))
          }
        }
      }
    }

    return(df)
  }

  # --- 加载差异基因结果 ---
  deg_results_from_file <- eventReactive(input$load_deg, {
    req(data_input$deg_file_data())

    showNotification("正在加载差异基因结果...", type = "message")

    df <- data_input$deg_file_data()
    cat("上传的文件列名:", paste(colnames(df), collapse = ", "), "\n")

    # 应用增强的列映射
    df <- enhanced_column_mapping(df)

    # 检查必要的列是否存在
    required_cols <- c("pvalue", "log2FoldChange")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
      showNotification(paste("缺少必要的列:", paste(missing_cols, collapse = ", ")), type = "error")
      showNotification("请确保上传的文件包含pvalue和log2FoldChange列，或使用以下列名之一:", type = "warning")
      showNotification("log2FoldChange: log2FC, avg_log2FC, logFC, log2_fold_change, log2fc, log2fc_adj", type = "message")
      showNotification("pvalue: p_val, p.value, P.Value, pvalue_adj", type = "message")
      return(NULL)
    }

    # 确保有padj列，如果没有则使用pvalue（但会标记警告）
    if (!"padj" %in% colnames(df)) {
      df$padj <- df$pvalue
      showNotification("⚠️ 警告：未找到校正p值（padj/FDR）列，将使用原始p值代替。", type = "warning")
      showNotification("建议：差异分析结果应包含多重假设检验校正后的p值。", type = "message")
      # 添加标记列，以便后续分析知道这是未校正的数据
      df$using_unadjusted_pval <- TRUE
    } else {
      df$using_unadjusted_pval <- FALSE
    }

    # 重命名列以匹配内部格式
    res <- df

    # 确保所有必要列都存在
    if (!"GeneID" %in% colnames(res)) {
      if ("SYMBOL" %in% colnames(res)) {
        res$GeneID <- res$SYMBOL
      } else {
        # 如果都没有，使用行名
        res$GeneID <- rownames(res)
      }
    }

    # 添加缺失的列
    if (!"baseMean" %in% colnames(res)) res$baseMean <- 1
    if (!"logCPM" %in% colnames(res)) res$logCPM <- 0

    # 确保SYMBOL列存在
    if (!"SYMBOL" %in% colnames(res)) {
      res$SYMBOL <- res$GeneID
    }

    # --- 差异状态判断 ---
    pval_col <- if(input$deg_pval_type == "p_val_adj") "padj" else "pvalue"

    res$Status <- ifelse(res[[pval_col]] < input$deg_pval_cutoff & abs(res$log2FoldChange) > input$deg_log2fc_cutoff,
                         ifelse(res$log2FoldChange > 0, "Up", "Down"), "Not DE")

    # TF 活性分析需要
    # 🔥 关键修复：确保t_stat不会产生Inf值
    res <- res %>%
      dplyr::mutate(
        # 限制pvalue的最小值，避免-log10(pvalue)过大
        pvalue_safe = pmax(pvalue, 1e-300),  # 防止log10(0) = Inf
        # 计算t_stat，并限制范围
        t_stat = -log10(pvalue_safe) * log2FoldChange
      ) %>%
      # 移除Inf和NA值
      dplyr::mutate(
        t_stat = ifelse(is.finite(t_stat), t_stat, NA)
      )

    cat(sprintf("📊 差异分析: %d 个基因的t_stat\n", sum(!is.na(res$t_stat))))
    cat(sprintf("📊 t_stat范围: %.2f 至 %.2f\n",
                min(res$t_stat, na.rm = TRUE),
                max(res$t_stat, na.rm = TRUE)))

    # --- 注释基因 ---
    anno <- data_input$annotate_genes(res$GeneID, input$deg_species)

    if (!is.null(anno)) {
      # 改进的GeneID清理逻辑 - 保留Ensembl ID的版本号
      clean_geneid <- res$GeneID
      clean_geneid <- trimws(clean_geneid)
      clean_geneid <- gsub("[\t\n\r]", "", clean_geneid)
      # 对于非Ensembl ID，才移除特殊字符
      non_ensembl <- !grepl("^ENS", clean_geneid, ignore.case = TRUE)
      clean_geneid[non_ensembl] <- gsub("[^[:alnum:]]", "", clean_geneid[non_ensembl])

      # 清理anno中的列名
      anno_clean <- anno
      if ("SYMBOL" %in% colnames(anno_clean)) {
        # 使用与原始数据相同的清理逻辑
        anno_clean$SYMBOL_CLEAN <- gsub("[^[:alnum:]]", "", anno_clean$SYMBOL)
        if (input$deg_species == "Mm") {
          anno_clean$SYMBOL_CLEAN <- sapply(anno_clean$SYMBOL_CLEAN, function(x) {
            if (grepl("^[A-Za-z]", x) && nchar(x) > 0) {
              paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
            } else {
              x
            }
          }, USE.NAMES = FALSE)
        } else {
          anno_clean$SYMBOL_CLEAN <- toupper(anno_clean$SYMBOL_CLEAN)
        }
      }

      # 初始化结果列
      if (!"SYMBOL" %in% colnames(res)) res$SYMBOL <- res$GeneID
      if (!"ENTREZID" %in% colnames(res)) res$ENTREZID <- NA

      # 🔥 修复：记录原始GeneID，用于后续显示
      res$Original_GeneID <- res$GeneID

      # 🔥 新的匹配策略：优先使用ENSEMBL列匹配Ensembl ID
      if ("ENSEMBL" %in% colnames(anno_clean)) {
        # 检查哪些基因是Ensembl ID
        is_ensembl_id <- grepl("^ENS", clean_geneid, ignore.case = TRUE)

        if (any(is_ensembl_id)) {
          # 直接通过ENSEMBL列匹配
          ensembl_match_idx <- match(clean_geneid[is_ensembl_id], anno_clean$ENSEMBL)
          matched_ensembl <- !is.na(ensembl_match_idx)

          if (any(matched_ensembl)) {
            ensembl_indices <- which(is_ensembl_id)
            indices_to_update <- ensembl_indices[matched_ensembl]
            res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[ensembl_match_idx[matched_ensembl]]
            res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[ensembl_match_idx[matched_ensembl]]
            cat("通过ENSEMBL列匹配成功:", sum(matched_ensembl), "个基因\n")
          }
        }
      }

      # 第一步：尝试SYMBOL匹配（对于非Ensembl ID或未匹配的基因）
      if ("SYMBOL" %in% colnames(anno_clean)) {
        # 只对尚未匹配ENTREZID的基因尝试SYMBOL匹配
        unmatched <- is.na(res$ENTREZID) | res$ENTREZID == ""

        if (any(unmatched)) {
          match_idx <- match(clean_geneid[unmatched], anno_clean$SYMBOL_CLEAN)
          matched_genes <- !is.na(match_idx)

          if (any(matched_genes)) {
            unmatched_indices <- which(unmatched)
            indices_to_update <- unmatched_indices[matched_genes]
            res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[match_idx[matched_genes]]
            res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[match_idx[matched_genes]]
            cat("通过SYMBOL列匹配成功:", sum(matched_genes), "个基因\n")
          }
        }

        # 🔥 第二步：对于仍未匹配的基因，尝试使用ENTREZID反向查询
        still_unmatched <- is.na(res$ENTREZID) | res$ENTREZID == ""
        if (any(still_unmatched)) {
          # 尝试通过ENTREZID匹配
          entrez_match_idx <- match(res$GeneID[still_unmatched], anno_clean$ENTREZID)
          matched_entrez <- !is.na(entrez_match_idx)

          if (any(matched_entrez)) {
            unmatched_indices <- which(still_unmatched)
            indices_to_update <- unmatched_indices[matched_entrez]
            res$SYMBOL[indices_to_update] <- anno_clean$SYMBOL[entrez_match_idx[matched_entrez]]
            res$ENTREZID[indices_to_update] <- anno_clean$ENTREZID[entrez_match_idx[matched_entrez]]
            cat("通过ENTREZID反向匹配成功:", sum(matched_entrez), "个基因\n")
          }
        }
      }
    } else {
      res$SYMBOL <- res$GeneID
      res$ENTREZID <- NA
      # 🔥 修复：记录原始GeneID
      res$Original_GeneID <- res$GeneID
    }

    if (!"SYMBOL" %in% colnames(res)) res$SYMBOL <- res$GeneID
    if (!"ENTREZID" %in% colnames(res)) res$ENTREZID <- NA

    # 🌟 过滤假基因
    res <- data_input$filter_pseudo_genes(res)

    # 改进的基因去重逻辑 - 保留统计显著性最高的基因
    if (any(duplicated(res$SYMBOL))) {
      # 按照p值和log2FoldChange的显著性排序
      res <- res %>%
        dplyr::arrange(SYMBOL, padj, abs(log2FoldChange)) %>%
        dplyr::distinct(SYMBOL, .keep_all = TRUE)

      # 记录去重信息
      n_duplicates_before <- nrow(res) - nrow(dplyr::distinct(res, SYMBOL))
      if (n_duplicates_before > 0) {
        cat(sprintf("移除了 %d 个重复的基因记录\n", n_duplicates_before))
      }
    }

    # 最终检查
    cat("最终数据列:", paste(colnames(res), collapse = ", "), "\n")
    cat("log2FoldChange类型:", class(res$log2FoldChange), "\n")
    cat("log2FoldChange范围:", range(res$log2FoldChange, na.rm=TRUE), "\n")

    return(res)
  })

  # 🆕 --- 加载芯片差异结果 ---
  chip_results_from_file <- eventReactive(input$load_chip, {
    req(data_input$chip_file_data())

    showNotification("正在加载芯片差异结果...", type = "message")

    df <- data_input$chip_file_data()
    cat("芯片文件列名:", paste(colnames(df), collapse = ", "), "\n")

    # 检查必要的列是否存在（芯片limma结果的列名）
    required_cols <- c("logFC", "P.Value", "SYMBOL", "ID")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
      showNotification(paste("缺少必要的列:", paste(missing_cols, collapse = ", ")), type = "error")
      showNotification("请确保上传的文件包含: logFC, AveExpr, t, P.Value, adj.P.Val, B, SYMBOL, ID", type = "warning")
      return(NULL)
    }

    # 转换为标准格式
    res <- data.frame(
      ID = as.character(df$ID),        # Entrez Gene ID (保持字符)
      SYMBOL = df$SYMBOL,              # 基因符号
      log2FoldChange = df$logFC,       # log2倍数变化
      pvalue = df$P.Value,             # 原始p值
      padj = df$adj.P.Val,             # BH校正p值
      baseMean = df$AveExpr,           # 平均表达
      t = df$t,                        # t统计量
      ENTREZID = as.numeric(as.character(df$ID)),  # 🔧 转换为numeric（clusterProfiler需要）
      GeneID = df$SYMBOL,              # 基因符号（用于兼容）
      Original_GeneID = df$SYMBOL,
      stringsAsFactors = FALSE
    )

    # 检查ENTREZID转换结果
    na_count <- sum(is.na(res$ENTREZID))
    if (na_count > 0) {
      cat(sprintf("⚠️ 警告: %d个基因的ENTREZID转换为NA（可能包含非数字ID）\n", na_count))
    }

    # --- 差异状态判断 ---
    pval_col <- if(input$chip_pval_type == "adj.P.Val") "padj" else "pvalue"

    res$Status <- ifelse(res[[pval_col]] < input$chip_pval_cutoff & abs(res$log2FoldChange) > input$chip_log2fc_cutoff,
                         ifelse(res$log2FoldChange > 0, "Up", "Down"), "Not DE")

    # TF 活性分析需要 - 计算t_stat
    res <- res %>%
      dplyr::mutate(
        pvalue_safe = pmax(pvalue, 1e-300),
        t_stat = -log10(pvalue_safe) * log2FoldChange
      ) %>%
      dplyr::mutate(
        t_stat = ifelse(is.finite(t_stat), t_stat, NA)
      )

    cat(sprintf("📊 芯片分析: %d 个基因的t_stat\n", sum(!is.na(res$t_stat))))

    # 芯片数据已经包含了ID和SYMBOL，不需要额外注释
    # 🔧 过滤掉ENTREZID为NA的基因（clusterProfiler需要有效的Entrez ID）
    before_filter <- nrow(res)
    res <- res[!is.na(res$ENTREZID), ]
    after_filter <- nrow(res)
    if (before_filter > after_filter) {
      cat(sprintf("⚠️ 过滤了 %d 个ENTREZID为NA的基因\n", before_filter - after_filter))
    }

    # 去重
    if (any(duplicated(res$SYMBOL))) {
      res <- res %>%
        dplyr::arrange(SYMBOL, padj, abs(log2FoldChange)) %>%
        dplyr::distinct(SYMBOL, .keep_all = TRUE)
    }

    cat("✅ 芯片数据加载完成:", nrow(res), "个基因\n")
    showNotification(sprintf("✅ 芯片数据加载完成: %d 个基因", nrow(res)), type = "message")

    return(res)
  })

  # --- 获取过滤后的表达矩阵基因列表（用于背景基因集） ---
  get_filtered_expr_genes <- reactive({
    if (input$data_source == "counts") {
      # 从原始数据开始
      req(data_input$raw_data(), input$control_group, input$treat_group)
      df <- data_input$raw_data()
      ctrl <- input$control_group
      trt <- input$treat_group
      df_use <- df[, c(ctrl, trt)]
      group <- factor(c(rep("Control", length(ctrl)), rep("Treatment", length(trt))),
                      levels = c("Control", "Treatment"))

      # 使用与deg_results()相同的过滤逻辑
      # 注意：filterByExpr()对两种方法使用相同的逻辑
      dge <- DGEList(counts = df_use, group = group)
      dge <- calcNormFactors(dge)
      keep <- filterByExpr(dge)
      filtered_genes <- rownames(dge)[keep]

      return(filtered_genes)
    } else {
      # 对于上传的差异基因文件，无法获取原始表达矩阵
      # 返回NULL，让富集分析模块使用默认背景
      return(NULL)
    }
  })

  # --- 统一的差异结果获取函数 ---
  get_deg_results <- reactive({
    if (input$data_source == "counts") {
      # 返回完整数据：deg_df + 表达矩阵 + 分组信息
      return(list(
        deg_df = deg_results(),
        background_genes = get_filtered_expr_genes(),
        expr_matrix = data_input$raw_data(),  # 完整表达矩阵
        ctrl_samples = input$ctrl_samples,     # 对照组样本
        trt_samples = input$trt_samples        # 处理组样本
      ))
    } else if (input$data_source == "deg") {
      # 上传差异文件时无法获取表达矩阵
      return(list(
        deg_df = deg_results_from_file(),
        background_genes = NULL,
        expr_matrix = NULL,
        ctrl_samples = NULL,
        trt_samples = NULL
      ))
    } else if (input$data_source == "chip") {
      # 🆕 芯片差异结果
      return(list(
        deg_df = chip_results_from_file(),
        background_genes = NULL,  # 芯片数据无法提供背景基因
        expr_matrix = NULL,
        ctrl_samples = NULL,
        trt_samples = NULL
      ))
    }
  })

  # --- 恢复主差异分析结果下载 ---
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("DEG_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(get_deg_results())
      write.csv(get_deg_results()$deg_df, file, row.names = FALSE)
    }
  )

  # --- 差异基因统计信息UI ---
  output$deg_summary <- renderUI({
    req(get_deg_results())
    data <- get_deg_results()$deg_df

    # 计算统计信息
    n_total <- nrow(data)
    n_up <- sum(data$Status == "Up", na.rm = TRUE)
    n_down <- sum(data$Status == "Down", na.rm = TRUE)
    n_not_de <- sum(data$Status == "Not DE", na.rm = TRUE)

    # 计算比例
    if (n_total > 0) {
      pct_up <- round(100 * n_up / n_total, 1)
      pct_down <- round(100 * n_down / n_total, 1)
      pct_not_de <- round(100 * n_not_de / n_total, 1)
    } else {
      pct_up <- pct_down <- pct_not_de <- 0
    }

    # 创建统计卡片
    tagList(
      tags$div(
        class = "row",
        style = "margin-bottom: 20px;",
        tags$div(
          class = "col-sm-3",
          tags$div(
            class = "card",
            style = "border-left: 4px solid #6c757d;",
            tags$div(class = "card-body", style = "padding: 15px;",
              tags$h6(class = "card-subtitle mb-2 text-muted", "总基因数"),
              tags$h3(class = "card-title mb-0", style = "color: #6c757d;",
                format(n_total, big.mark = ",")
              ),
              tags$small(class = "text-muted", paste0("100%"))
            )
          )
        ),
        tags$div(
          class = "col-sm-3",
          tags$div(
            class = "card",
            style = "border-left: 4px solid #28a745;",
            tags$div(class = "card-body", style = "padding: 15px;",
              tags$h6(class = "card-subtitle mb-2 text-muted", "上调基因"),
              tags$h3(class = "card-title mb-0", style = "color: #28a745;",
                format(n_up, big.mark = ",")
              ),
              tags$small(class = "text-muted", paste0(pct_up, "%"))
            )
          )
        ),
        tags$div(
          class = "col-sm-3",
          tags$div(
            class = "card",
            style = "border-left: 4px solid #dc3545;",
            tags$div(class = "card-body", style = "padding: 15px;",
              tags$h6(class = "card-subtitle mb-2 text-muted", "下调基因"),
              tags$h3(class = "card-title mb-0", style = "color: #dc3545;",
                format(n_down, big.mark = ",")
              ),
              tags$small(class = "text-muted", paste0(pct_down, "%"))
            )
          )
        ),
        tags$div(
          class = "col-sm-3",
          tags$div(
            class = "card",
            style = "border-left: 4px solid #17a2b8;",
            tags$div(class = "card-body", style = "padding: 15px;",
              tags$h6(class = "card-subtitle mb-2 text-muted", "非显著"),
              tags$h3(class = "card-title mb-0", style = "color: #17a2b8;",
                format(n_not_de, big.mark = ",")
              ),
              tags$small(class = "text-muted", paste0(pct_not_de, "%"))
            )
          )
        )
      )
    )
  })

  output$deg_table <- DT::renderDataTable({
    req(get_deg_results())
    data_to_display <- get_deg_results()$deg_df

    # 只格式化存在的列
    numeric_cols <- c("log2FoldChange", "pvalue", "padj", "t_stat")
    existing_numeric_cols <- numeric_cols[numeric_cols %in% colnames(data_to_display)]

    if (length(existing_numeric_cols) > 0) {
      DT::datatable(data_to_display, options = list(scrollX=T, pageLength=10), rownames=F) %>%
        formatRound(existing_numeric_cols, 4)
    } else {
      DT::datatable(data_to_display, options = list(scrollX=T, pageLength=10), rownames=F)
    }
  })

  # --- 自定义基因显示 ---
  custom_genes <- reactiveVal(NULL)

  observeEvent(input$show_custom_genes, {
    req(input$custom_genes_input)

    # 解析用户输入的基因
    genes <- strsplit(input$custom_genes_input, ",")[[1]]
    genes <- trimws(genes)  # 去除空格
    genes <- genes[genes != ""]  # 去除空字符串

    if (length(genes) > 0) {
      custom_genes(genes)
      showNotification(paste("已设置显示", length(genes), "个自定义基因"), type = "message")
    } else {
      custom_genes(NULL)
      showNotification("请输入有效的基因名称", type = "warning")
    }
  })

  # 清除自定义基因
  observeEvent(input$clear_custom_genes, {
    custom_genes(NULL)
    updateTextInput(session, "custom_genes_input", value = "")
    showNotification("已清除自定义基因", type = "message")
  })


  # --- 火山图 ---
  output$interactive_volcano <- renderPlotly({
    req(get_deg_results())
    res_data <- get_deg_results()
    res <- res_data$deg_df  # 获取实际的数据框

    # 添加调试信息
    cat("火山图数据检查:\n")
    cat("数据类型:", class(res), "\n")
    cat("数据列名:", paste(colnames(res), collapse = ", "), "\n")
    if ("log2FoldChange" %in% colnames(res)) {
      cat("log2FoldChange类型:", class(res$log2FoldChange), "\n")
    }

    # 根据数据来源选择p值类型
    if (input$data_source == "counts") {
      pval_col <- input$pval_type
    } else {
      pval_col <- if(input$deg_pval_type == "p_val_adj") "padj" else "pvalue"
    }

    # 使用用户选择的Y轴类型
    y_axis_col <- input$y_axis_type

    # 检查log2FoldChange列
    if (!("log2FoldChange" %in% colnames(res) && is.numeric(res$log2FoldChange))) {
      showNotification("错误：log2FoldChange列不存在或不是数值类型", type = "error")
      showNotification(paste("当前列名:", paste(colnames(res), collapse = ", ")), type = "message")
      return(NULL)
    }

    # 安全计算-log10值，处理非数值和NA值
    if (y_axis_col %in% colnames(res) && is.numeric(res[[y_axis_col]])) {
      # 确保数值有效且大于0（log10需要正数）
      valid_values <- res[[y_axis_col]]

      # 使用机器最小正值代替0，避免log10(0)的问题
      min_positive <- .Machine$double.xmin  # 约为2.2e-308
      valid_values[valid_values <= 0 & !is.na(valid_values)] <- min_positive
      valid_values[is.na(valid_values)] <- NA

      res$y_value <- -log10(valid_values)

      # 检查是否有有效的y值
      if (all(is.na(res$y_value))) {
        showNotification(paste("错误：所有", y_axis_col, "值无效（<=0或NA），无法绘制火山图"), type = "error")
        return(NULL)
      }

      # 如果有极小值被替换，给出警告
      n_replaced <- sum(res[[y_axis_col]] <= 0 & !is.na(res[[y_axis_col]]))
      if (n_replaced > 0) {
        showNotification(sprintf("注意：有 %d 个p值为0或负值的基因被替换为最小正值", n_replaced), type = "message")
      }
    } else {
      showNotification(paste("错误：列", y_axis_col, "不存在或不是数值类型"), type = "error")
      return(NULL)
    }

    color_map <- c("Not DE"="#95a5a6", "Up"=input$up_color, "Down"=input$down_color)

    txt_col <- if(input$theme_toggle) "#00e0ff" else "black"

    # 创建基础火山图
    p <- plot_ly(res,
                 x = ~log2FoldChange, y = ~y_value, color = ~Status,
                 colors = color_map,
                 text = ~SYMBOL, type = 'scatter', mode = 'markers',
                 marker = list(size = input$point_size, opacity = input$point_alpha),
                 hoverinfo = 'text',
                 hovertext = ~paste("Gene:", SYMBOL,
                                   "<br>log2FC:", round(log2FoldChange, 3),
                                   "<br>-log10(", y_axis_col, "):", round(y_value, 3),
                                   "<br>Status:", Status)) %>%
      layout(
        xaxis = list(
          title = "log2(Fold Change)",
          range = c(input$x_axis_min, input$x_axis_max),
          titlefont = list(size = input$axis_title_size),
          tickfont = list(size = input$axis_label_size),
          showgrid = input$show_grid,
          gridcolor = if(input$show_grid) "#ddd" else "transparent",
          gridwidth = 1
        ),
        yaxis = list(
          title = paste0("-log10(", y_axis_col, ")"),
          titlefont = list(size = input$axis_title_size),
          tickfont = list(size = input$axis_label_size),
          showgrid = input$show_grid,
          gridcolor = if(input$show_grid) "#ddd" else "transparent",
          gridwidth = 1
        ),
        font = list(color = txt_col),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)"
      )

    # 添加自定义基因标签
    if (!is.null(custom_genes())) {
      selected_genes <- custom_genes()

      # 在结果中查找这些基因
      gene_data <- res[res$SYMBOL %in% selected_genes | res$GeneID %in% selected_genes, ]

      if (nrow(gene_data) > 0) {
        # 添加基因标签
        p <- p %>%
          add_annotations(
            x = gene_data$log2FoldChange,
            y = gene_data$y_value,
            text = gene_data$SYMBOL,
            xref = "x",
            yref = "y",
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = 1,
            arrowwidth = 1,
            arrowcolor = input$gene_label_color,
            ax = 20,
            ay = -40,
            font = list(
              size = input$gene_label_size,
              color = input$gene_label_color,
              family = "Arial",
              weight = if(input$gene_label_bold) "bold" else "normal"
            ),
            bgcolor = "rgba(255,255,255,0.8)",
            bordercolor = input$gene_label_color,
            borderwidth = 1,
            borderpad = 4,
            opacity = 0.8
          )
      }
    }

    p
  })

  # --- 静态火山图用于导出 ---
  volcano_static_plot <- reactive({
    req(get_deg_results())
    res_data <- get_deg_results()
    res <- res_data$deg_df  # 获取实际的数据框

    # 根据数据来源选择p值类型
    if (input$data_source == "counts") {
      pval_col <- input$pval_type
    } else {
      pval_col <- if(input$deg_pval_type == "p_val_adj") "padj" else "pvalue"
    }

    # 使用用户选择的Y轴类型
    y_axis_col <- input$y_axis_type

    # 检查log2FoldChange列
    if (!("log2FoldChange" %in% colnames(res) && is.numeric(res$log2FoldChange))) {
      showNotification("错误：log2FoldChange列不存在或不是数值类型", type = "error")
      showNotification(paste("当前列名:", paste(colnames(res), collapse = ", ")), type = "message")
      return(NULL)
    }

    # 安全计算-log10值，处理非数值和NA值
    if (y_axis_col %in% colnames(res) && is.numeric(res[[y_axis_col]])) {
      # 确保数值有效且大于0（log10需要正数）
      valid_values <- res[[y_axis_col]]

      # 使用机器最小正值代替0，避免log10(0)的问题
      min_positive <- .Machine$double.xmin  # 约为2.2e-308
      valid_values[valid_values <= 0 & !is.na(valid_values)] <- min_positive
      valid_values[is.na(valid_values)] <- NA

      res$y_value <- -log10(valid_values)

      # 检查是否有有效的y值
      if (all(is.na(res$y_value))) {
        showNotification(paste("错误：所有", y_axis_col, "值无效（<=0或NA），无法绘制火山图"), type = "error")
        return(NULL)
      }

      # 如果有极小值被替换，给出警告
      n_replaced <- sum(res[[y_axis_col]] <= 0 & !is.na(res[[y_axis_col]]))
      if (n_replaced > 0) {
        showNotification(sprintf("注意：有 %d 个p值为0或负值的基因被替换为最小正值", n_replaced), type = "message")
      }
    } else {
      showNotification(paste("错误：列", y_axis_col, "不存在或不是数值类型"), type = "error")
      return(NULL)
    }

    # 设置颜色
    res$color <- ifelse(res$Status == "Up", input$up_color,
                       ifelse(res$Status == "Down", input$down_color, "#95a5a6"))

    # 创建ggplot火山图（与交互图保持一致的大小比例）
    p <- ggplot(res, aes(x = log2FoldChange, y = y_value, color = Status)) +
      geom_point(alpha = input$point_alpha, size = input$point_size) +
      scale_color_manual(values = c("Up" = input$up_color, "Down" = input$down_color, "Not DE" = "#95a5a6")) +
      labs(
        x = "log2(Fold Change)",
        y = paste0("-log10(", y_axis_col, ")"),
        title = "Volcano Plot"
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = input$axis_title_size),
        axis.text = element_text(size = input$axis_label_size),
        legend.title = element_text(size = input$axis_title_size),
        legend.text = element_text(size = input$axis_label_size),
        plot.title = element_text(size = input$axis_title_size + 2, hjust = 0.5),
        panel.grid.major = element_line(
          color = if(input$show_grid) "gray" else "transparent",
          linewidth = if(input$show_grid) 0.5 else 0
        ),
        panel.grid.minor = element_line(
          color = if(input$show_grid) "gray" else "transparent",
          linewidth = if(input$show_grid) 0.25 else 0
        )
      ) +
      xlim(input$x_axis_min, input$x_axis_max)

    # 添加自定义基因标签（与交互图保持一致的大小）
    if (!is.null(custom_genes())) {
      selected_genes <- custom_genes()
      gene_data <- res[res$SYMBOL %in% selected_genes | res$GeneID %in% selected_genes, ]

      if (nrow(gene_data) > 0) {
        p <- p +
          geom_text_repel(
            data = gene_data,
            aes(label = SYMBOL),
            size = input$gene_label_size,
            color = input$gene_label_color,
            fontface = if(input$gene_label_bold) "bold" else "plain",
            box.padding = 0.5,
            point.padding = 0.3,
            max.overlaps = Inf
          )
      }
    }

    return(p)
  })

  # --- 火山图导出 ---
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", Sys.Date(), ".", input$export_format)
    },
    content = function(file) {
      req(volcano_static_plot())

      if (input$export_format == "png") {
        png(file, width = input$export_width, height = input$export_height, units = "in", res = 300)
      } else if (input$export_format == "pdf") {
        pdf(file, width = input$export_width, height = input$export_height)
      } else if (input$export_format == "svg") {
        svg(file, width = input$export_width, height = input$export_height)
      }

      print(volcano_static_plot())
      dev.off()
    }
  )

  # 返回差异分析结果
  return(get_deg_results)
}