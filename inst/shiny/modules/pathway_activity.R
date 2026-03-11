# =====================================================
# 通路活性分析模块
# =====================================================

#' 通路活性分析服务器模块
#' @param input Shiny输入
#' @param output Shiny输出
#' @param session Shiny会话
#' @param deg_results 差异分析结果
#' @param kegg_results KEGG富集分析结果
pathway_activity_server <- function(input, output, session, deg_results, kegg_results) {

  # =====================================================
  # 1. 运行通路活性分析
  # =====================================================

  pathway_activity_results <- eventReactive(input$run_pathway_activity, {
    req(deg_results(), kegg_results())

    # 获取选择的方法
    method <- input$pathway_method

    showNotification(paste("正在推断通路活性 (方法:", toupper(method), ")..."),
                    type = "message")

    # 从差异分析结果中提取数据
    deg_data <- deg_results()
    deg_res <- deg_data$deg_df

    # 检查是否有数据
    if (is.null(deg_res) || nrow(deg_res) == 0) {
      showNotification("通路活性分析失败: 差异分析结果不存在", type = "error")
      return(NULL)
    }

    # 1. 构建表达矩阵
    # 根据选择的方法决定使用完整表达矩阵还是 log2FoldChange

    # 🆕 检查是否有完整表达矩阵（用于 AUCell/GSVA）
    has_expr_matrix <- !is.null(deg_data$expr_matrix)

    # 🔍 方法验证
    if (method %in% c("aucell", "gsva")) {
      if (!has_expr_matrix) {
        showNotification(
          "⚠️ AUCell/GSVA 方法需要从 counts 数据运行差异分析。",
          "当前使用的是上传的差异基因文件，无法获得完整表达矩阵。",
          "自动切换到 ULM 方法。",
          type = "warning",
          duration = 10
        )
        method <- "ulm"
      } else {
        n_samples <- ncol(deg_data$expr_matrix)
        cat(sprintf("📊 使用完整表达矩阵: %d 基因 × %d 样本\n",
                    nrow(deg_data$expr_matrix), n_samples))

        if (n_samples < 3) {
          showNotification(
            paste("💡 提示: AUCell/GSVA 建议至少有 3 个样本。当前样本数:", n_samples),
            type = "message",
            duration = 5
          )
        }
      }
    }

    # 使用 log2FoldChange 作为基因表达变化度量 (ULM/WMEAN)
    # 检查是否有ENTREZID列（用于与KEGG匹配）
    has_entrez <- "ENTREZID" %in% colnames(deg_res)

    # 根据KEGG使用的ID类型，提前准备匹配用的数据框
    # 注意：这里需要保留ENTREZID列，因为如果KEGG使用ENTREZID，我们之后需要它
    if (has_entrez) {
      stats_df <- deg_res %>%
        filter(!is.na(SYMBOL), !is.na(log2FoldChange), !is.na(ENTREZID)) %>%
        group_by(SYMBOL) %>%
        filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
        ungroup() %>%
        distinct(SYMBOL, .keep_all = TRUE) %>%
        select(SYMBOL, ENTREZID, log2FoldChange) %>%  # 保留ENTREZID以备后用
        mutate(ENTREZID = as.character(ENTREZID)) %>%
        filter(is.finite(log2FoldChange))
    } else {
      stats_df <- deg_res %>%
        filter(!is.na(SYMBOL), !is.na(log2FoldChange)) %>%
        group_by(SYMBOL) %>%
        filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
        ungroup() %>%
        distinct(SYMBOL, .keep_all = TRUE) %>%
        select(SYMBOL, log2FoldChange) %>%
        filter(is.finite(log2FoldChange))
    }

    if (nrow(stats_df) < 5) {
      showNotification(
        paste0("通路分析失败: 有效基因数量 (", nrow(stats_df), ") 不足"),
        type = "error",
        duration = 10
      )
      return(NULL)
    }

    # 2. 构建 KEGG 通路网络
    # 直接从 kegg_results reactive 函数获取
    kegg_data <- tryCatch({
      kegg_results()
    }, error = function(e) {
      showNotification("通路活性分析失败: 无法获取 KEGG 富集结果", type = "error")
      return(NULL)
    })

    if (is.null(kegg_data) || nrow(kegg_data) == 0) {
      showNotification("通路活性分析失败: 请先运行 KEGG 富集分析", type = "error", duration = 10)
      return(NULL)
    }

    # 📊 调试信息：检查KEGG数据格式
    cat("\n=== KEGG数据调试 ===\n")
    cat(sprintf("KEGG结果行数: %d\n", nrow(kegg_data)))
    cat(sprintf("KEGG基因ID示例 (前5个): %s\n",
                paste(head(kegg_data$geneID, 5), collapse = ", ")))

    # 解析一个geneID看格式
    sample_geneid <- kegg_data$geneID[1]
    parsed_genes <- strsplit(sample_geneid, "/")[[1]]
    cat(sprintf("示例geneID: %s\n", sample_geneid))
    cat(sprintf("解析后的基因: %s\n", paste(head(parsed_genes, 5), collapse = ", ")))
    cat(sprintf("基因类型: %s\n", class(parsed_genes[1])))

    # 🔍 关键判断：KEGG使用的是SYMBOL还是ENTREZID？
    # 从示例基因判断：如果是数字开头→ENTREZID，如果是字母开头→SYMBOL
    first_gene <- parsed_genes[1]
    kegg_uses_symbol <- grepl("^[A-Za-z]", first_gene)

    cat(sprintf("KEGG使用: %s\n", ifelse(kegg_uses_symbol, "SYMBOL", "ENTREZID")))

    # 解析基因ID (从 geneID 列)
    # 根据KEGG使用的ID类型选择匹配方式
    if (kegg_uses_symbol) {
      # KEGG使用SYMBOL，我们也用SYMBOL匹配
      expr_df <- stats_df %>%
        select(SYMBOL, log2FC = log2FoldChange)
      rownames_label <- "SYMBOL"
      # 保留log2FoldChange列作为唯一的数据列
      mat_input <- stats_df %>%
        select(SYMBOL, log2FoldChange) %>%  # 只选择需要的列
        column_to_rownames(var = "SYMBOL") %>%
        as.matrix()
      # 确保列名是 log2FoldChange
      if (ncol(mat_input) == 1) {
        colnames(mat_input) <- "log2FoldChange"
      }
    } else {
      # KEGG使用ENTREZID
      if (has_entrez) {
        expr_df <- stats_df %>%
          select(ENTREZID, log2FC = log2FoldChange)
        rownames_label <- "ENTREZID"
        mat_input <- stats_df %>%
          select(ENTREZID, log2FoldChange) %>%  # 只选择需要的列
          column_to_rownames(var = "ENTREZID") %>%
          as.matrix()
        # 确保列名是 log2FoldChange
        if (ncol(mat_input) == 1) {
          colnames(mat_input) <- "log2FoldChange"
        }
      } else {
        # 没有ENTREZID，只能用SYMBOL
        expr_df <- stats_df %>%
          select(SYMBOL, log2FC = log2FoldChange)
        rownames_label <- "SYMBOL"
        mat_input <- stats_df %>%
          select(SYMBOL, log2FoldChange) %>%  # 只选择需要的列
          column_to_rownames(var = "SYMBOL") %>%
          as.matrix()
        # 确保列名是 log2FoldChange
        if (ncol(mat_input) == 1) {
          colnames(mat_input) <- "log2FoldChange"
        }
      }
    }

    # 📊 调试信息：检查expr_df格式
    cat(sprintf("\n=== 表达数据调试 ===\n"))
    cat(sprintf("表达矩阵行数: %d\n", nrow(expr_df)))
    cat(sprintf("使用的ID类型: %s\n", rownames_label))
    cat(sprintf("ID示例 (前5个): %s\n",
                paste(head(expr_df[[rownames_label]], 5), collapse = ", ")))
    cat(sprintf("ID类型: %s\n", class(expr_df[[rownames_label]][1])))

    # 构建通路网络
    # MOR（mode of regulation）应为「通路对基因的生物学调控方向」：
    # 若未知（KEGG 仅提供通路-基因隶属关系），应全部设为 1，由基因的 log2FC 决定正负贡献。
    # 错误做法：mor = ifelse(log2FC>0, 1, -1) 会使 贡献=log2FC*mor 恒为正，导致所有通路得分>0。
    pathway_net <- kegg_data %>%
      filter(!is.na(geneID)) %>%
      mutate(
        pathway = Description,
        genes = strsplit(geneID, "/")
      ) %>%
      unnest(genes) %>%
      left_join(expr_df, by = c("genes" = rownames_label)) %>%
      filter(!is.na(log2FC)) %>%
      group_by(pathway) %>%
      mutate(mor = 1L) %>%   # 无方向信息时统一为 1，通路得分 = 基因 log2FC 的聚合，可正可负
      ungroup() %>%
      select(
        source = pathway,
        target = genes,
        mor
      )

    # 📊 诊断信息：检查匹配情况
    mat_rownames <- rownames(mat_input)
    unique_targets <- unique(pathway_net$target)
    matched_targets <- sum(unique_targets %in% mat_rownames)

    cat(sprintf("📊 通路网络初步构建: %d 通路, %d 相互关系\n",
                length(unique(pathway_net$source)), nrow(pathway_net)))
    cat(sprintf("📊 表达矩阵中的基因数: %d\n", length(mat_rownames)))
    cat(sprintf("📊 通路网络中的唯一基因数: %d\n", length(unique_targets)))
    cat(sprintf("📊 匹配的基因数: %d (%.1f%%)\n",
                matched_targets, 100*matched_targets/length(unique_targets)))

    # 最后的过滤：确保target在表达矩阵中
    pathway_net <- pathway_net %>%
      filter(target %in% mat_rownames) %>%
      group_by(source) %>%
      filter(n() >= input$pathway_minsize) %>%
      ungroup()

    if (nrow(pathway_net) == 0) {
      showNotification("通路活性分析失败: 没有找到足够的通路-基因映射", type = "error", duration = 10)
      cat("❌ 错误: 过滤后通路网络为空\n")
      cat("可能原因:\n")
      cat("  1. ENTREZID类型不匹配 (字符 vs 数字)\n")
      cat("  2. 最小基因集大小设置过大\n")
      cat("  3. KEGG通路基因与差异基因重叠太少\n")
      return(NULL)
    }

    # 📊 诊断信息
    mor_distribution <- table(pathway_net$mor)
    cat(sprintf("📊 通路网络构建完成: %d 通路, %d 相互关系\n",
                length(unique(pathway_net$source)), nrow(pathway_net)))
    cat(sprintf("   MOR分布: 激活=%d, 抑制=%d\n",
                sum(pathway_net$mor == 1, na.rm = TRUE),
                sum(pathway_net$mor == -1, na.rm = TRUE)))
    cat(sprintf("📊 表达矩阵维度: %d 基因 x %d 样本\n",
                nrow(mat_input), ncol(mat_input)))

    # 检查log2FoldChange分布
    log2fc_values <- as.numeric(mat_input[, "log2FoldChange"])
    log2fc_range <- range(log2fc_values, na.rm = TRUE)
    cat(sprintf("📊 log2FoldChange范围: [%.3f, %.3f]\n",
                log2fc_range[1], log2fc_range[2]))
    cat(sprintf("📊 正值基因数: %d, 负值基因数: %d\n",
                sum(log2fc_values > 0),
                sum(log2fc_values < 0)))

    # 3. 运行通路活性推断
    tryCatch({
      # 🆕 根据方法选择合适的数据源和算法
      result_df <- switch(method,
        "ulm" = {
          decoupleR::run_ulm(
            mat = mat_input,
            net = pathway_net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$pathway_minsize
          )
        },
        "wmean" = {
          decoupleR::run_wmean(
            mat = mat_input,
            net = pathway_net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$pathway_minsize
          )
        },
        "aucell" = {
          # 🆕 AUCell: 使用完整表达矩阵
          if (!has_expr_matrix) {
            # 这不会执行，因为前面已经检查并降级了
            stop("AUCell requires expression matrix")
          }

          # 准备表达矩阵（确保行名是基因）
          expr_mat <- deg_data$expr_matrix
          if (!is.null(rownames(expr_mat))) {
            rownames(expr_mat) <- expr_mat[[1]]
            expr_mat <- expr_mat[, -1, drop=FALSE]
          }

          decoupleR::run_aucell(
            mat = expr_mat,
            net = pathway_net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$pathway_minsize
          )
        },
        "gsva" = {
          # 🆕 GSVA: 使用完整表达矩阵
          if (!has_expr_matrix) {
            # 这不会执行，因为前面已经检查并降级了
            stop("GSVA requires expression matrix")
          }

          # 准备表达矩阵（确保行名是基因）
          expr_mat <- deg_data$expr_matrix
          if (!is.null(rownames(expr_mat))) {
            rownames(expr_mat) <- expr_mat[[1]]
            expr_mat <- expr_mat[, -1, drop=FALSE]
          }

          decoupleR::run_gsva(
            mat = expr_mat,
            net = pathway_net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$pathway_minsize
          )
        },
        {
          # 默认使用 ULM
          decoupleR::run_ulm(
            mat = mat_input,
            net = pathway_net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$pathway_minsize
          )
        }
      )

      # 处理结果
      if (is.data.frame(result_df)) {
        df <- result_df
      } else if (is.list(result_df)) {
        if ("statistic" %in% names(result_df)) {
          df <- result_df$statistic
        } else {
          df <- as.data.frame(result_df)
        }
      } else {
        df <- as.data.frame(result_df)
      }

      # 添加排名
      df <- df %>%
        mutate(rnk = NA)

      msk_pos <- df$score > 0
      df[msk_pos, 'rnk'] <- rank(-df[msk_pos, 'score'])

      msk_neg <- df$score < 0
      df[msk_neg, 'rnk'] <- rank(-abs(df[msk_neg, 'score']))

      # 添加方法信息
      df$method <- toupper(method)

      # 📊 结果诊断信息
      score_summary <- summary(df$score)
      cat(sprintf("\n📊 通路活性推断结果统计:\n"))
      cat(sprintf("   总通路数: %d\n", nrow(df)))
      cat(sprintf("   活跃通路 (score>0): %d (%.1f%%)\n",
                  sum(df$score > 0), 100*mean(df$score > 0)))
      cat(sprintf("   抑制通路 (score<0): %d (%.1f%%)\n",
                  sum(df$score < 0), 100*mean(df$score < 0)))
      cat(sprintf("   Score范围: [%.4f, %.4f]\n",
                  min(df$score), max(df$score)))
      cat(sprintf("   Score中位数: %.4f\n", median(df$score)))

      showNotification(paste("通路活性推断完成! (方法:", toupper(method),
                               ", 通路数:", nrow(df), ")"),
                      type = "message")

      return(df)

    }, error = function(e) {
      showNotification(paste("通路活性分析失败 (", toupper(method), "):", e$message),
                      type = "error")
      return(NULL)
    })
  })

  # =====================================================
  # 3. 可视化：柱状图
  # =====================================================

  output$pathway_activity_bar_plot <- renderPlot({
    req(pathway_activity_results())

    df_acts <- pathway_activity_results()

    n_pathways <- input$pathway_top_n

    pathways_to_plot <- df_acts %>%
      arrange(rnk) %>%
      head(n_pathways) %>%
      pull(source)

    f_pathway_acts <- df_acts %>%
      filter(source %in% pathways_to_plot)

    txt_col <- if(input$theme_toggle) "white" else "black"
    grid_col <- if(input$theme_toggle) "#444444" else "#cccccc"

    # 获取用户设置的字体大小
    user_fontsize <- input$pathway_bar_fontsize

    ggplot(f_pathway_acts, aes(x = reorder(source, score), y = score)) +
      geom_bar(aes(fill = score), stat = "identity") +
      scale_fill_gradient2(
        low = input$pathway_inactive_col,
        high = input$pathway_active_col,
        mid = "whitesmoke",
        midpoint = 0,
        name = "活性分数"
      ) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = txt_col) +
      coord_flip() +
      theme_minimal() +
      labs(
        x = "KEGG 通路",
        y = paste("活性分数 (", toupper(df_acts$method[1]), " Score", ")"),
        title = paste("Top", n_pathways, "KEGG 通路活性变化 (T vs C)")
      ) +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5, size = user_fontsize + 4),
        axis.title = element_text(color = txt_col, face = "bold", size = user_fontsize + 2),
        axis.text.x = element_text(size = user_fontsize, face = "bold", color = txt_col),
        axis.text.y = element_text(size = user_fontsize, face = "bold", color = txt_col),
        legend.text = element_text(color = txt_col, size = user_fontsize),
        legend.title = element_text(color = txt_col, size = user_fontsize + 1, face = "bold"),
        axis.line = element_line(color = txt_col),
        panel.grid.major = element_line(color = grid_col),
        panel.grid.minor = element_line(color = grid_col)
      )
  })

  # =====================================================
  # 4. 可视化：热图
  # =====================================================

  output$pathway_activity_heatmap <- renderPlot({
    req(pathway_activity_results())

    df_acts <- pathway_activity_results()

    n_pathways <- input$pathway_top_n

    # 选择 top 通路
    top_pathways <- df_acts %>%
      arrange(rnk) %>%
      head(n_pathways) %>%
      pull(source)

    heatmap_data <- df_acts %>%
      filter(source %in% top_pathways) %>%
      select(source, score) %>%
      spread(key = source, value = score)

    if (ncol(heatmap_data) < 2) {
      return(NULL)
    }

    # 转置以正确显示
    heatmap_matrix <- as.matrix(heatmap_data)
    rownames(heatmap_matrix) <- "Activity"

    txt_col <- if(input$theme_toggle) "white" else "black"

    # 获取用户设置的字体大小
    user_fontsize <- input$pathway_heatmap_fontsize

    pheatmap(
      heatmap_matrix,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      color = colorRampPalette(c(input$pathway_inactive_col, "white",
                                  input$pathway_active_col))(100),
      main = paste("Top", n_pathways, "通路活性热图"),
      fontsize = user_fontsize,           # 使用用户设置的字体大小
      fontsize_row = user_fontsize,        # 行名（Activity）
      fontsize_col = user_fontsize - 2     # 列名（通路名称）稍小一点
    )
  })

  # =====================================================
  # 5. 结果表格
  # =====================================================

  output$pathway_activity_table <- DT::renderDataTable({
    req(pathway_activity_results())

    df <- pathway_activity_results() %>%
      select(source, score, p_value, rnk) %>%
      rename(Pathway = source, Score = score, P.Value = p_value, Rank = rnk) %>%
      arrange(Rank)

    DT::datatable(df, selection = 'single',
                 options = list(scrollX=T, pageLength=15),
                 rownames=F) %>%
      formatRound(c("Score", "P.Value"), 4)
  })

  # =====================================================
  # 6. 下载结果
  # =====================================================

  output$download_pathway_results <- downloadHandler(
    filename = function() {
      paste0("Pathway_Activity_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pathway_activity_results())
      df <- pathway_activity_results() %>%
        select(source, score, p_value, rnk) %>%
        rename(Pathway = source, Score = score, P.Value = p_value, Rank = rnk)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # 6.1 下载热图 PNG
  output$download_pathway_heatmap_png <- downloadHandler(
    filename = function() {
      paste0("Pathway_Activity_Heatmap_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(pathway_activity_results())

      df_acts <- pathway_activity_results()
      n_pathways <- input$pathway_top_n

      # 选择 top 通路
      top_pathways <- df_acts %>%
        arrange(rnk) %>%
        head(n_pathways) %>%
        pull(source)

      heatmap_data <- df_acts %>%
        filter(source %in% top_pathways) %>%
        select(source, score) %>%
        spread(key = source, value = score)

      if (ncol(heatmap_data) < 2) {
        return(NULL)
      }

      # 转置以正确显示
      heatmap_matrix <- as.matrix(heatmap_data)
      rownames(heatmap_matrix) <- "Activity"

      # 获取用户设置的字体大小
      user_fontsize <- input$pathway_heatmap_fontsize

      # 保存为 PNG
      png(file, width = 1200, height = 800, res = 150)

      pheatmap(
        heatmap_matrix,
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        display_numbers = TRUE,
        number_format = "%.2f",
        color = colorRampPalette(c(input$pathway_inactive_col, "white",
                                    input$pathway_active_col))(100),
        main = paste("Top", n_pathways, "通路活性热图"),
        fontsize = user_fontsize,
        fontsize_row = user_fontsize,
        fontsize_col = user_fontsize - 2
      )

      dev.off()
    }
  )

  # 6.2 下载热图 SVG
  output$download_pathway_heatmap_svg <- downloadHandler(
    filename = function() {
      paste0("Pathway_Activity_Heatmap_", Sys.Date(), ".svg")
    },
    content = function(file) {
      req(pathway_activity_results())

      df_acts <- pathway_activity_results()
      n_pathways <- input$pathway_top_n

      # 选择 top 通路
      top_pathways <- df_acts %>%
        arrange(rnk) %>%
        head(n_pathways) %>%
        pull(source)

      heatmap_data <- df_acts %>%
        filter(source %in% top_pathways) %>%
        select(source, score) %>%
        spread(key = source, value = score)

      if (ncol(heatmap_data) < 2) {
        return(NULL)
      }

      # 转置以正确显示
      heatmap_matrix <- as.matrix(heatmap_data)
      rownames(heatmap_matrix) <- "Activity"

      # 获取用户设置的字体大小
      user_fontsize <- input$pathway_heatmap_fontsize

      # 保存为 SVG (使用 svglight 设备)
      # 检查是否安装了 svglight 包，如果没有则使用备用方案
      if (requireNamespace("svglight", quietly = TRUE)) {
        svglight::svglight(file, width = 12, height = 8)

        pheatmap(
          heatmap_matrix,
          cluster_rows = FALSE,
          cluster_cols = TRUE,
          display_numbers = TRUE,
          number_format = "%.2f",
          color = colorRampPalette(c(input$pathway_inactive_col, "white",
                                      input$pathway_active_col))(100),
          main = paste("Top", n_pathways, "通路活性热图"),
          fontsize = user_fontsize,
          fontsize_row = user_fontsize,
          fontsize_col = user_fontsize - 2
        )

        dev.off()
      } else {
        # 如果没有 svglight，使用 recordPlot 并保存为 RDS
        # 创建临时 PDF，然后提示用户转换
        temp_pdf <- tempfile(fileext = ".pdf")
        pdf(temp_pdf, width = 12, height = 8)

        pheatmap(
          heatmap_matrix,
          cluster_rows = FALSE,
          cluster_cols = TRUE,
          display_numbers = TRUE,
          number_format = "%.2f",
          color = colorRampPalette(c(input$pathway_inactive_col, "white",
                                      input$pathway_active_col))(100),
          main = paste("Top", n_pathways, "通路活性热图"),
          fontsize = user_fontsize,
          fontsize_row = user_fontsize,
          fontsize_col = user_fontsize - 2
        )

        dev.off()

        # 复制 PDF 到目标文件（用户可以手动转换为 SVG）
        file.copy(temp_pdf, file, overwrite = TRUE)
        unlink(temp_pdf)

        showNotification("提示：已保存为 PDF 格式。如需 SVG，请安装 svglight 包：install.packages('svglight')",
                        type = "warning", duration = 10)
      }
    }
  )

  # =====================================================
  # 6.3 下载柱状图 PNG
  # =====================================================

  output$download_pathway_bar_png <- downloadHandler(
    filename = function() {
      paste0("Pathway_Activity_Barplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(pathway_activity_results())

      df_acts <- pathway_activity_results()
      n_pathways <- input$pathway_top_n

      pathways_to_plot <- df_acts %>%
        arrange(rnk) %>%
        head(n_pathways) %>%
        pull(source)

      f_pathway_acts <- df_acts %>%
        filter(source %in% pathways_to_plot)

      txt_col <- "black"
      grid_col <- "#cccccc"

      # 获取用户设置的字体大小
      user_fontsize <- input$pathway_bar_fontsize

      # 创建 ggplot 对象
      p <- ggplot(f_pathway_acts, aes(x = reorder(source, score), y = score)) +
        geom_bar(aes(fill = score), stat = "identity") +
        scale_fill_gradient2(
          low = input$pathway_inactive_col,
          high = input$pathway_active_col,
          mid = "whitesmoke",
          midpoint = 0,
          name = "活性分数"
        ) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = txt_col) +
        coord_flip() +
        theme_minimal() +
        labs(
          x = "KEGG 通路",
          y = paste("活性分数 (", toupper(df_acts$method[1]), " Score", ")"),
          title = paste("Top", n_pathways, "KEGG 通路活性变化 (T vs C)")
        ) +
        theme(
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5, size = user_fontsize + 4),
          axis.title = element_text(color = txt_col, face = "bold", size = user_fontsize + 2),
          axis.text.x = element_text(size = user_fontsize, face = "bold", color = txt_col),
          axis.text.y = element_text(size = user_fontsize, face = "bold", color = txt_col),
          legend.text = element_text(color = txt_col, size = user_fontsize),
          legend.title = element_text(color = txt_col, size = user_fontsize + 1, face = "bold"),
          axis.line = element_line(color = txt_col),
          panel.grid.major = element_line(color = grid_col),
          panel.grid.minor = element_line(color = grid_col)
        )

      # 保存为 PNG
      ggsave(
        filename = file,
        plot = p,
        device = "png",
        width = 12,
        height = 8,
        dpi = 300,
        bg = "white"
      )
    }
  )

  # 6.4 下载柱状图 SVG
  output$download_pathway_bar_svg <- downloadHandler(
    filename = function() {
      paste0("Pathway_Activity_Barplot_", Sys.Date(), ".svg")
    },
    content = function(file) {
      req(pathway_activity_results())

      df_acts <- pathway_activity_results()
      n_pathways <- input$pathway_top_n

      pathways_to_plot <- df_acts %>%
        arrange(rnk) %>%
        head(n_pathways) %>%
        pull(source)

      f_pathway_acts <- df_acts %>%
        filter(source %in% pathways_to_plot)

      txt_col <- "black"
      grid_col <- "#cccccc"

      # 获取用户设置的字体大小
      user_fontsize <- input$pathway_bar_fontsize

      # 创建 ggplot 对象
      p <- ggplot(f_pathway_acts, aes(x = reorder(source, score), y = score)) +
        geom_bar(aes(fill = score), stat = "identity") +
        scale_fill_gradient2(
          low = input$pathway_inactive_col,
          high = input$pathway_active_col,
          mid = "whitesmoke",
          midpoint = 0,
          name = "活性分数"
        ) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = txt_col) +
        coord_flip() +
        theme_minimal() +
        labs(
          x = "KEGG 通路",
          y = paste("活性分数 (", toupper(df_acts$method[1]), " Score", ")"),
          title = paste("Top", n_pathways, "KEGG 通路活性变化 (T vs C)")
        ) +
        theme(
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5, size = user_fontsize + 4),
          axis.title = element_text(color = txt_col, face = "bold", size = user_fontsize + 2),
          axis.text.x = element_text(size = user_fontsize, face = "bold", color = txt_col),
          axis.text.y = element_text(size = user_fontsize, face = "bold", color = txt_col),
          legend.text = element_text(color = txt_col, size = user_fontsize),
          legend.title = element_text(color = txt_col, size = user_fontsize + 1, face = "bold"),
          axis.line = element_line(color = txt_col),
          panel.grid.major = element_line(color = grid_col),
          panel.grid.minor = element_line(color = grid_col)
        )

      # 保存为 SVG - 使用 R 内置的 svg() 函数
      svg(file, width = 12, height = 8)
      print(p)
      dev.off()
    }
  )

  # =====================================================
  # 7. KEGG 模块联动说明
  # =====================================================
  # 此模块直接接收 kegg_results 作为参数
  # 当 KEGG 富集分析完成后，通路活性分析会自动获取结果
  # 用户只需切换到"🛤️ 通路活性"标签即可开始分析

  # =====================================================
  # 8. 方法说明
  # =====================================================

  output$pathway_method_info <- renderUI({
    method <- input$pathway_method

    info_text <- switch(method,
      "aucell" = list(
        title = "AUCell 方法",
        desc = "基于基因集富集分析，计算通路在每个样本中的富集分数",
        features = c("• 适用于单细胞和 bulk RNA-seq",
                    "• 计算基因集的 AUC 分数",
                    "• 无需预设分组",
                    "• 对噪声数据较稳健")
      ),
      "gsva" = list(
        title = "GSVA 方法",
        desc = "基因集变异分析，提供样本层面的通路活性评分",
        features = c("• 样本层面的连续评分",
                    "• 无需预设分组",
                    "• 适用于时间序列分析",
                    "• 可识别渐进变化")
      ),
      "ulm" = list(
        title = "ULM 方法 (单变量线性模型)",
        desc = "基于线性模型的通路活性推断",
        features = c("• 快速且稳定",
                    "• 适合大规模数据",
                    "• 考虑基因调控模式",
                    "• 提供统计显著性评估")
      ),
      "wmean" = list(
        title = "WMEAN 方法 (加权平均)",
        desc = "基于加权平均的通路活性评分",
        features = c("• 计算简单快速",
                    "• 权重基于调控模式",
                    "• 适合探索性分析",
                    "• 结果易于解释")
      )
    )

    tags$div(
      class = "alert alert-info",
      style = "margin: 10px 0;",
      h5(info_text$title, style = "color: #0c5460; margin-top: 0;"),
      p(info_text$desc, style = "margin: 5px 0;"),
      tags$ul(style = "margin: 5px 0; padding-left: 20px;",
        lapply(info_text$features, function(f) {
          tags$li(f, style = "margin: 3px 0;")
        })
      )
    )
  })

  # =====================================================
  # 9. 统计摘要
  # =====================================================

  output$pathway_summary <- renderText({
    req(pathway_activity_results())

    df <- pathway_activity_results()

    n_total <- nrow(df)
    n_active <- sum(df$score > 0, na.rm = TRUE)
    n_inactive <- sum(df$score < 0, na.rm = TRUE)
    n_significant <- sum(df$p_value < 0.05, na.rm = TRUE)

    paste0(
      "📊 通路活性统计 | ",
      "总通路数: ", n_total, " | ",
      "激活: ", n_active, " (", round(100*n_active/n_total, 1), "%) | ",
      "抑制: ", n_inactive, " (", round(100*n_inactive/n_total, 1), "%) | ",
      "显著 (P<0.05): ", n_significant
    )
  })

  # =====================================================
  # 10. 算法说明模块 UI 输出
  # =====================================================

  # 10.1 什么是通路活性分析
  output$pathway_algorithm_intro <- renderUI({
    tags$div(
      class = "algorithm-section",
      tags$p(
        style = "font-size: 14px; line-height: 1.6;",
        "通路活性分析是一种基于",
        tags$strong("基因表达变化", style = "color: #e74c3c;"),
        "和",
        tags$strong("通路-基因映射", style = "color: #3498db;"),
        "的系统生物学方法，用于推断特定生物通路在实验条件下的活跃程度。"
      ),
      tags$div(
        style = "background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("🎯 核心思想", style = "color: #2c3e50;"),
        tags$ul(
          tags$li("如果某个通路被", tags$span("激活", style = "color: #e74c3c; font-weight: bold;"), "，该通路中的上调基因应该整体表现出较高的表达变化"),
          tags$li("如果某个通路被", tags$span("抑制", style = "color: #3498db; font-weight: bold;"), "，该通路中的下调基因应该整体表现出较低的表达变化"),
          tags$li("通过统计模型量化基因表达变化与通路活性之间的关联")
        )
      ),
      tags$h5("✨ 主要用途", style = "color: #2c3e50; margin-top: 15px;"),
      tags$ul(
        tags$li("识别实验中", tags$strong("显著激活或抑制的信号通路")),
        tags$li("理解疾病发生发展的", tags$strong("分子机制")),
        tags$li("发现潜在的", tags$strong("药物靶点")),
        tags$li("验证", tags$strong("KEGG 富集分析结果", style = "color: #27ae60;"), "的生物学意义")
      )
    )
  })

  # 10.2 ULM 方法原理详解
  output$pathway_ulm_explanation <- renderUI({
    tags$div(
      class = "algorithm-section",
      tags$p(
        style = "font-size: 14px; line-height: 1.6;",
        tags$strong("ULM (Univariate Linear Model)", style = "font-size: 16px;"),
        "即单变量线性模型，是通过线性回归分析推断通路活性的统计方法。"
      ),

      # 核心公式
      tags$div(
        style = "background: #e8f4f8; padding: 20px; border-radius: 8px; margin: 15px 0; border-left: 5px solid #3498db;",
        tags$h5("📐 核心公式", style = "color: #2c3e50;"),
        tags$p(
          tags$code(
            style = "font-size: 18px; background: white; padding: 10px; border-radius: 4px; display: inline-block;",
            "Y = β₀ + β₁ × MOR + ε"
          )
        ),
        tags$ul(
          tags$li(tags$code("Y"), ": 基因的", tags$strong("log2FoldChange"), "（表达变化）"),
          tags$li(tags$code("MOR"), ": 调控模式（", tags$span("+1", style = "color: #e74c3c; font-weight: bold;"), " = 上调，",
                  tags$span("-1", style = "color: #3498db; font-weight: bold;"), " = 下调）"),
          tags$li(tags$code("β₁"), ": 回归系数 = ", tags$strong("通路活性分数", style = "color: #9b59b6;")),
          tags$li(tags$code("ε"), ": 残差（随机误差）")
        )
      ),

      # 工作流程
      tags$div(
        style = "background: #fef5e7; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("⚙️ 计算流程", style = "color: #2c3e50;"),
        tags$ol(
          tags$li("从差异分析结果提取基因的", tags$strong("log2FoldChange")),
          tags$li("根据 log2FC 正负计算每个基因的", tags$strong("MOR", "（+1 或 -1）")),
          tags$li("对于每个通路，提取其包含的所有基因"),
          tags$li("拟合线性模型：", tags$code("log2FC ~ MOR")),
          tags$li("提取回归系数", tags$code("β₁"), "作为该通路的", tags$strong("活性分数"))
        )
      ),

      # 结果解读
      tags$div(
        style = "background: #eafaf1; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("📊 结果解读", style = "color: #2c3e50;"),
        tags$table(
          style = "width: 100%; border-collapse: collapse;",
          tags$thead(
            tags$tr(
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #27ae60; color: white;", "Score (β₁)"),
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #27ae60; color: white;", "P < 0.05"),
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #27ae60; color: white;", "结论")
            )
          ),
          tags$tbody(
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center; color: #e74c3c; font-weight: bold;", "> 0"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center;", "✓"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$span("通路显著激活", style = "color: #e74c3c; font-weight: bold;"))
            ),
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center; color: #3498db; font-weight: bold;", "< 0"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center;", "✓"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$span("通路显著抑制", style = "color: #3498db; font-weight: bold;"))
            ),
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center;", "≈ 0"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px; text-align: center;", "✗"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "通路不活跃或变化不显著")
            )
          )
        )
      ),

      # 优势
      tags$div(
        style = "background: #f4ecf7; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("✅ 方法优势", style = "color: #2c3e50;"),
        tags$ul(
          tags$li(tags$strong("简单直观"), ": 基于经典线性回归，易于理解"),
          tags$li(tags$strong("统计严谨"), ": 提供 p 值评估显著性"),
          tags$li(tags$strong("计算快速"), ": 几秒完成分析，适合大规模数据"),
          tags$li(tags$strong("适用场景"), ": 差异表达分析结果（log2FoldChange）")
        )
      )
    )
  })

  # 10.3 结果解读指南
  output$pathway_result_guide <- renderUI({
    tags$div(
      class = "algorithm-section",
      tags$p("运行通路活性分析后，您将看到以下结果：", style = "font-size: 14px;"),

      # 控制台输出示例
      tags$div(
        style = "background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; font-family: 'Courier New', monospace; font-size: 12px;",
        tags$h5("💻 控制台输出示例", style = "color: #3498db; margin-top: 0;"),
        tags$pre(
          "📊 通路活性推断结果统计:\n",
          "   总通路数: 275\n",
          "   活跃通路 (score>0): 137 (49.8%)\n",
          "   抑制通路 (score<0): 138 (50.2%)\n",
          "   Score范围: [-0.2345, 0.2678]\n",
          "   Score中位数: 0.0012"
        )
      ),

      # 柱状图解读
      tags$div(
        style = "background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("📊 柱状图解读", style = "color: #2c3e50;"),
        tags$ul(
          tags$li(tags$span(style = "display: inline-block; width: 15px; height: 15px; background: #e74c3c; margin-right: 5px;"),
                  "红色柱子（Score > 0）→ 激活的通路"),
          tags$li(tags$span(style = "display: inline-block; width: 15px; height: 15px; background: #3498db; margin-right: 5px;"),
                  "蓝色柱子（Score < 0）→ 抑制的通路"),
          tags$li("柱子", tags$strong("越长"), "表示活性", tags$strong("越显著")),
          tags$li("按", tags$code("Rank"), "排序，排名越靠前说明活性越显著")
        )
      ),

      # 数据表解读
      tags$div(
        style = "background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("📋 数据表解读", style = "color: #2c3e50;"),
        tags$table(
          style = "width: 100%; border-collapse: collapse; font-size: 13px;",
          tags$thead(
            tags$tr(
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "列名"),
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "含义"),
              tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "如何使用")
            )
          ),
          tags$tbody(
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("Pathway")),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "KEGG 通路名称"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "识别感兴趣的通路")
            ),
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("Score")),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "活性分数（β₁）"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "> 0 激活，< 0 抑制，绝对值越大越显著")
            ),
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("P.Value")),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "统计显著性"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "< 0.05 认为显著")
            ),
            tags$tr(
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("Rank")),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "活性排名"),
              tags$td(style = "border: 1px solid #ddd; padding: 8px;", "数字越小，活性越显著")
            )
          )
        )
      ),

      # 关键检查点
      tags$div(
        style = "background: #fff3cd; padding: 15px; border-radius: 5px; margin: 10px 0; border-left: 5px solid #f39c12;",
        tags$h5("⚠️ 结果质量检查", style = "color: #2c3e50;"),
        tags$ul(
          tags$li(tags$strong("MOR 分布"), ": 激活≈50%, 抑制≈50% 为理想",
                  tags$br(), tags$small("如果激活=100%, 抑制=0% 说明基因匹配有问题")),
          tags$li(tags$strong("Score 范围"), ": 通常在 [-0.5, +0.5]",
                  tags$br(), tags$small("如果 >1.0 或 <-1.0 需要检查数据")),
          tags$li(tags$strong("显著通路比例"), ": P < 0.05 的通路应占 20-50%",
                  tags$br(), tags$small("如果 >80% 或 <5% 需要检查数据质量"))
        )
      )
    )
  })

  # 10.4 实际应用案例
  output$pathway_example_use_case <- renderUI({
    tags$div(
      class = "algorithm-section",

      tags$div(
        style = "background: #e8f6f3; padding: 20px; border-radius: 8px; margin: 15px 0; border-left: 5px solid #16a085;",
        tags$h5("🔬 案例: 炎症反应研究", style = "color: #2c3e50;"),
        tags$p(style = "font-size: 14px; margin-bottom: 10px;",
               "研究问题：", tags$strong("某药物对小鼠炎症模型的影响")),

        tags$h6("分析流程", style = "color: #16a085;"),
        tags$ol(
          tags$li("上传 RNA-seq 数据（药物处理组 vs 对照组）"),
          tags$li("运行差异分析，得到 2,345 个差异表达基因"),
          tags$li("运行 KEGG 富集分析，识别出 275 个富集通路"),
          tags$li("运行通路活性分析（ULM 方法）")
        ),

        tags$h6("关键发现", style = "color: #16a085; margin-top: 10px;"),
        tags$div(
          style = "background: white; padding: 10px; border-radius: 4px; margin: 5px 0;",
          tags$table(
            style = "width: 100%; font-size: 13px;",
            tags$thead(
              tags$tr(
                tags$th(style = "padding: 5px; text-align: left;", "通路"),
                tags$th(style = "padding: 5px; text-align: center;", "Score"),
                tags$th(style = "padding: 5px; text-align: center;", "P值"),
                tags$th(style = "padding: 5px; text-align: left;", "解读")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$td(style = "padding: 5px;", "TNF signaling pathway"),
                tags$td(style = "padding: 5px; text-align: center; color: #e74c3c; font-weight: bold;", "+0.1567"),
                tags$td(style = "padding: 5px; text-align: center;", "0.0001"),
                tags$td(style = "padding: 5px;", tags$span("显著激活", style = "color: #e74c3c; font-weight: bold;"))
              ),
              tags$tr(
                tags$td(style = "padding: 5px;", "NF-kappa B signaling"),
                tags$td(style = "padding: 5px; text-align: center; color: #e74c3c; font-weight: bold;", "+0.1234"),
                tags$td(style = "padding: 5px; text-align: center;", "0.0003"),
                tags$td(style = "padding: 5px;", tags$span("显著激活", style = "color: #e74c3c; font-weight: bold;"))
              ),
              tags$tr(
                tags$td(style = "padding: 5px;", "Oxidative phosphorylation"),
                tags$td(style = "padding: 5px; text-align: center; color: #3498db; font-weight: bold;", "-0.2345"),
                tags$td(style = "padding: 5px; text-align: center;", "0.0023"),
                tags$td(style = "padding: 5px;", tags$span("显著抑制", style = "color: #3498db; font-weight: bold;"))
              )
            )
          )
        ),

        tags$h6("生物学结论", style = "color: #16a085; margin-top: 10px;"),
        tags$p("药物显著激活了炎症相关通路（TNF、NF-κB），抑制了氧化磷酸化，说明",
               tags$strong("药物引发炎症反应", style = "color: #e74c3c;"), "并影响能量代谢。")
      ),

      # 应用价值
      tags$div(
        style = "background: #fef5e7; padding: 15px; border-radius: 5px; margin: 10px 0;",
        tags$h5("💡 通路活性分析的价值", style = "color: #2c3e50;"),
        tags$ul(
          tags$li(tags$strong("超越富集分析"), ": 不仅知道通路富集，还知道通路是激活还是抑制"),
          tags$li(tags$strong("量化活性程度"), ": Score 值可以比较不同通路的活性强度"),
          tags$li(tags$strong("统计显著性"), ": P 值帮助排除随机噪声"),
          tags$li(tags$strong("指导后续实验"), ": 选择显著激活/抑制的通路进行验证")
        )
      )
    )
  })

  # 10.5 常见问题
  output$pathway_faq <- renderUI({
    tags$div(
      class = "algorithm-section",

      tags$div(
        style = "margin-bottom: 20px;",
        tags$h5("❓ Q1: AUCell/GSVA 和 ULM/WMEAN 有什么区别？", style = "color: #e74c3c;"),
        tags$p(tags$strong("A:"), "它们使用不同的数据源和算法：",
              tags$br(), tags$br(),
              tags$table(
                style = "width: 100%; border-collapse: collapse;",
                tags$thead(
                  tags$tr(
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "方法"),
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "数据源"),
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "适用场景")
                  )
                ),
                tags$tbody(
                  tags$tr(
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$strong("AUCell/GSVA")),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "完整表达矩阵（多样本）"),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "从 counts 数据运行时可用")
                  ),
                  tags$tr(
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$strong("ULM/WMEAN")),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "log2FoldChange（单列）"),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "所有情况都可用")
                  )
                )
              ),
              tags$br(),
              tags$span("💡 建议：", style = "color: #27ae60; font-weight: bold;"),
              "如果有完整表达矩阵（从 counts 运行），推荐使用 AUCell/GSVA 获得更准确的结果。")
      ),

      tags$div(
        style = "margin-bottom: 20px;",
        tags$h5("❓ Q2: Score 的正常范围是多少？", style = "color: #e74c3c;"),
        tags$p(tags$strong("A:"), "理论上 `(-∞, +∞)`，但实际通常在", tags$code("[-0.5, +0.5]"),
              "范围内。",
              tags$br(), tags$br(),
              tags$ul(
                tags$li(tags$strong("Score > 0.3"), ": 强激活"),
                tags$li(tags$strong("Score < -0.3"), ": 强抑制"),
                tags$li(tags$strong("-0.1 < Score < 0.1"), ": 弱活性")
              ),
              tags$small("注意：Score 绝对值大小受数据质量影响，需结合 p-value 判断。", style = "color: #7f8c8d;"))
      ),

      tags$div(
        style = "margin-bottom: 20px;",
        tags$h5("❓ Q3: 如何调整参数获得更好的结果？", style = "color: #e74c3c;"),
        tags$p(tags$strong("A:"), "关键是调整", tags$code("最小基因集大小"),
              tags$br(), tags$br(),
              tags$table(
                style = "width: 100%; border-collapse: collapse;",
                tags$tr(
                  tags$td(style = "border: 1px solid #ddd; padding: 8px; background: #ecf0f1;", tags$strong("参数")),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px; background: #ecf0f1;", tags$strong("效果")),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px; background: #ecf0f1;", tags$strong("推荐值"))
                ),
                tags$tr(
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("minsize = 3")),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", "更敏感，检测更多通路"),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", "探索性分析")
                ),
                tags$tr(
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("minsize = 5-10")),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", "平衡敏感性和可靠性"),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$span("✅ 默认推荐", style = "color: #27ae60;"))
                ),
                tags$tr(
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$code("minsize = 15+")),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", "更保守，减少噪声"),
                  tags$td(style = "border: 1px solid #ddd; padding: 8px;", "验证性分析")
                )
              ))
      ),

      tags$div(
        style = "margin-bottom: 20px;",
        tags$h5("❓ Q4: 通路活性和 KEGG 富集分析有什么区别？", style = "color: #e74c3c;"),
        tags$p(tags$strong("A:"), "两者互补但不同：",
              tags$br(), tags$br(),
              tags$table(
                style = "width: 100%; border-collapse: collapse;",
                tags$thead(
                  tags$tr(
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "方面"),
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "KEGG 富集分析"),
                    tags$th(style = "border: 1px solid #ddd; padding: 8px; background: #34495e; color: white;", "通路活性分析")
                  )
                ),
                tags$tbody(
                  tags$tr(
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$strong("目的")),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "识别差异基因富集在哪些通路"),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "判断通路是激活还是抑制")
                  ),
                  tags$tr(
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$strong("输出")),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "P值、富集分数"),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "Score（活性分数）、P值")
                  ),
                  tags$tr(
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", tags$strong("解读")),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "通路是否显著富集"),
                    tags$td(style = "border: 1px solid #ddd; padding: 8px;", "通路激活（+）或抑制（-）")
                  )
                )
              ),
              tags$br(),
              tags$span("💡 建议：", style = "color: #27ae60; font-weight: bold;"),
              "先运行 KEGG 富集分析识别重要通路，再用通路活性分析判断它们的激活状态。")
      ),

      tags$div(
        style = "margin-bottom: 20px;",
        tags$h5("❓ Q5: 为什么所有通路都是激活（100%）？", style = "color: #e74c3c;"),
        tags$p(tags$strong("A:"), "这通常表示基因匹配有问题。检查控制台输出：",
              tags$br(), tags$br(),
              tags$code("MOR分布: 激活=5834, 抑制=0"),
              tags$br(), tags$br(),
              "可能原因：",
              tags$ul(
                tags$li("KEGG 使用的 ID 类型（SYMBOL vs ENTREZID）与差异分析结果不匹配"),
                tags$li("log2FoldChange 计算有问题"),
                tags$li("基因去重逻辑有误")
              ),
              tags$span("✅ 解决方案：", style = "color: #27ae60;"),
      "检查控制台调试信息，确认 KEGG 和 DEG 结果的 ID 类型一致。")
      )
    )
  })

  # 返回通路活性结果供其他模块使用
  return(pathway_activity_results)
}
