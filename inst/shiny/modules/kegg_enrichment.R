# =====================================================
# KEGG富集分析模块
# =====================================================

kegg_enrichment_server <- function(input, output, session, deg_results) {

  # =====================================================
  # 辅助函数：清理基因符号
  # =====================================================
  clean_gene_symbols <- function(gene_symbols, species_code) {
    # 清理基因符号：去除空格、特殊字符，标准化大小写
    cleaned <- trimws(gene_symbols)  # 去除首尾空格
    cleaned <- gsub("[\t\n\r]", "", cleaned)  # 去除空白字符

    # 去除版本号（如.1, .2等）
    cleaned <- gsub("\\.[0-9]+$", "", cleaned)

    # 去除常见的假基因后缀
    cleaned <- gsub("-ps$", "", cleaned, ignore.case = TRUE)
    cleaned <- gsub("-rs$", "", cleaned, ignore.case = TRUE)
    cleaned <- gsub("-as$", "", cleaned, ignore.case = TRUE)

    # 识别并处理ENSEMBL ID
    # ENSEMBL ID模式：ENS(MUS)?G[0-9]+（人类：ENSG，小鼠：ENSMUSG）
    is_ensembl_id <- grepl("^ENS(MUS)?G[0-9]+$", cleaned, ignore.case = TRUE)

    # 根据物种和ID类型标准化大小写
    if (species_code == "mmu") {
      # 小鼠基因处理
      cleaned <- sapply(seq_along(cleaned), function(i) {
        gene <- cleaned[i]

        if (is_ensembl_id[i]) {
          # ENSEMBL ID：全部大写
          return(toupper(gene))
        } else if (grepl("^[A-Za-z]", gene)) {
          # 基因符号：首字母大写，其余小写
          return(paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene)))))
        } else {
          # 其他情况（如数字ID）
          return(gene)
        }
      }, USE.NAMES = FALSE)
    } else {
      # 人类基因：全部大写（包括ENSEMBL ID和基因符号）
      cleaned <- toupper(cleaned)
    }

    # 去除连字符、点等特殊字符（保留字母和数字）
    # 注意：对于ENSEMBL ID，这可能会去除有效的版本号，但我们已经处理了版本号
    cleaned <- gsub("[^[:alnum:]]", "", cleaned)

    return(cleaned)
  }

  # =====================================================
  # 辅助函数：识别基因ID类型
  # =====================================================
  identify_gene_id_types <- function(gene_ids, species_code) {
    # 初始化结果列表
    result <- list(
      ensembl_ids = character(0),
      gene_symbols = character(0),
      entrez_ids = character(0),
      other_ids = character(0)
    )

    for (gene in gene_ids) {
      # 检查是否是ENSEMBL ID
      if (grepl("^ENS(MUS)?G[0-9]+$", gene, ignore.case = TRUE)) {
        result$ensembl_ids <- c(result$ensembl_ids, gene)
      }
      # 检查是否是ENTREZID（纯数字）
      else if (grepl("^[0-9]+$", gene)) {
        result$entrez_ids <- c(result$entrez_ids, gene)
      }
      # 检查是否是基因符号（以字母开头）
      else if (grepl("^[A-Za-z]", gene)) {
        result$gene_symbols <- c(result$gene_symbols, gene)
      }
      # 其他类型
      else {
        result$other_ids <- c(result$other_ids, gene)
      }
    }

    return(result)
  }

  # =====================================================
  # 辅助函数：智能基因符号转换
  # =====================================================
  smart_gene_conversion <- function(gene_ids, db_obj, target_column = "ENTREZID") {
    # 尝试不同的keytype来转换基因ID
    keytypes_to_try <- c("SYMBOL", "ALIAS", "ENSEMBL", "ENTREZID")

    # 记录调试信息
    debug_info <- list()
    debug_info$input_count <- length(gene_ids)
    debug_info$input_samples <- head(gene_ids, 10)
    debug_info$attempts <- list()

    for (keytype in keytypes_to_try) {
      tryCatch({
        # 先检查哪些基因ID在当前keytype中有效
        valid_keys <- keys(db_obj, keytype = keytype)
        matched_ids <- gene_ids[gene_ids %in% valid_keys]

        # 记录尝试信息
        attempt_info <- list(
          keytype = keytype,
          valid_keys_count = length(valid_keys),
          matched_count = length(matched_ids),
          matched_samples = head(matched_ids, 5)
        )
        debug_info$attempts[[keytype]] <- attempt_info

        if (length(matched_ids) > 0) {
          # 尝试转换匹配的基因ID
          converted <- AnnotationDbi::mapIds(
            db_obj,
            keys = matched_ids,
            column = target_column,
            keytype = keytype,
            multiVals = "first"
          )

          # 返回成功转换的结果
          successful <- converted[!is.na(converted)]
          if (length(successful) > 0) {
            debug_info$final_keytype <- keytype
            debug_info$success_count <- length(successful)
            debug_info$failed_count <- length(matched_ids) - length(successful)

            # 打印调试信息（在开发环境中）
            if (Sys.getenv("SHINY_DEBUG") == "TRUE") {
              cat("\n=== 基因转换调试信息 ===\n")
              cat("输入基因数量:", debug_info$input_count, "\n")
              cat("输入基因示例:", paste(debug_info$input_samples, collapse=", "), "\n")
              cat("成功使用的keytype:", keytype, "\n")
              cat("匹配的基因数量:", length(matched_ids), "\n")
              cat("成功转换的基因数量:", length(successful), "\n")
              cat("失败的基因数量:", debug_info$failed_count, "\n")
            }

            return(list(
              converted = successful,
              keytype_used = keytype,
              matched_count = length(matched_ids),
              success_count = length(successful),
              debug_info = debug_info
            ))
          }
        }
      }, error = function(e) {
        # 记录错误信息
        debug_info$attempts[[keytype]]$error <- e$message
        # 继续尝试下一个keytype
        NULL
      })
    }

    # 如果所有keytype都失败，返回详细的调试信息
    debug_info$all_failed <- TRUE

    # 打印详细的调试信息
    if (Sys.getenv("SHINY_DEBUG") == "TRUE") {
      cat("\n=== 基因转换失败调试信息 ===\n")
      cat("输入基因数量:", debug_info$input_count, "\n")
      cat("输入基因示例:", paste(debug_info$input_samples, collapse=", "), "\n")
      cat("\n尝试的keytype结果:\n")
      for (keytype in keytypes_to_try) {
        if (!is.null(debug_info$attempts[[keytype]])) {
          attempt <- debug_info$attempts[[keytype]]
          cat("  ", keytype, ":\n")
          cat("    有效key数量:", attempt$valid_keys_count, "\n")
          cat("    匹配数量:", attempt$matched_count, "\n")
          if (!is.null(attempt$matched_samples)) {
            cat("    匹配示例:", paste(attempt$matched_samples, collapse=", "), "\n")
          }
          if (!is.null(attempt$error)) {
            cat("    错误:", attempt$error, "\n")
          }
        }
      }
    }

    return(list(
      converted = NULL,
      keytype_used = NULL,
      matched_count = 0,
      success_count = 0,
      debug_info = debug_info,
      error_message = "所有keytype尝试都失败了"
    ))
  }

  # =====================================================
  # KEGG 模块
  # =====================================================
  kegg_data_processed <- eventReactive(input$run_kegg, {
    req(deg_results())

    # 从deg_results中提取差异分析结果和背景基因
    deg_data <- deg_results()
    res_df <- deg_data$deg_df
    background_genes <- deg_data$background_genes

    target_status <- switch(input$kegg_direction, "Up" = "Up", "Down" = "Down", "All" = c("Up", "Down"))

    # 清理基因符号
    if (!is.null(background_genes) && length(background_genes) > 0) {
      background_genes <- clean_gene_symbols(background_genes, input$kegg_species)
    }

    # 获取ENTREZID
    ids <- res_df %>% dplyr::filter(Status %in% target_status & !is.na(ENTREZID)) %>% dplyr::pull(ENTREZID)

    if(length(ids) == 0) {
      showNotification("无有效ENTREZID，请检查基因注释结果", type="error")
      return(NULL)
    }

    tryCatch({
      # 准备背景基因集（如果可用）
      universe <- NULL
      if(!is.null(background_genes) && length(background_genes) > 0) {
        # 将背景基因符号转换为ENTREZID（使用智能转换）
        db_pkg <- if(input$kegg_species == "mmu") "org.Mm.eg.db" else "org.Hs.eg.db"
        if(require(db_pkg, character.only = TRUE)) {
          db_obj <- get(db_pkg)

          # 先清理基因符号
          cleaned_background <- clean_gene_symbols(background_genes, input$kegg_species)

          # 使用智能转换函数
          conversion_result <- smart_gene_conversion(cleaned_background, db_obj, "ENTREZID")

          # 检查转换结果
          if(!is.null(conversion_result$converted) && length(conversion_result$converted) > 0) {
            bg_entrez <- conversion_result$converted
            universe <- bg_entrez

            # 显示成功信息
            if(!is.null(conversion_result$keytype_used)) {
              showNotification(paste("使用", length(universe), "个检测到的基因作为KEGG分析背景基因集（通过",
                                    conversion_result$keytype_used, "转换）"), type = "message")
            } else {
              showNotification(paste("使用", length(universe), "个检测到的基因作为KEGG分析背景基因集"), type = "message")
            }

            # 如果有失败的情况，显示警告
            if(conversion_result$matched_count > conversion_result$success_count) {
              failed_count <- conversion_result$matched_count - conversion_result$success_count
              showNotification(paste("警告：", failed_count, "个背景基因无法转换为ENTREZID"), type = "warning")
            }
          } else {
            # 转换失败，显示详细的错误信息
            error_msg <- "背景基因转换失败"

            if(!is.null(conversion_result$error_message)) {
              error_msg <- paste(error_msg, ": ", conversion_result$error_message)
            }

            # 提供具体的建议
            if(length(cleaned_background) > 0) {
              sample_genes <- head(cleaned_background, 5)
              error_msg <- paste0(error_msg, "\n示例基因：", paste(sample_genes, collapse=", "))

              # 分析基因ID类型并提供具体建议
              id_types <- identify_gene_id_types(sample_genes, input$kegg_species)
              error_msg <- paste0(error_msg, "\n\n检测到的ID类型分析：")

              if(length(id_types$ensembl_ids) > 0) {
                error_msg <- paste0(error_msg, "\n• ENSEMBL ID: ", length(id_types$ensembl_ids), "个")
                error_msg <- paste0(error_msg, "\n  示例：", paste(head(id_types$ensembl_ids, 3), collapse=", "))
                error_msg <- paste0(error_msg, "\n  建议：这些是ENSEMBL ID，不是基因符号。")
                error_msg <- paste0(error_msg, "\n  请使用基因符号（如Trp53）或确保数据库包含这些ENSEMBL ID")
              }

              if(length(id_types$gene_symbols) > 0) {
                error_msg <- paste0(error_msg, "\n• 基因符号: ", length(id_types$gene_symbols), "个")
                error_msg <- paste0(error_msg, "\n  示例：", paste(head(id_types$gene_symbols, 3), collapse=", "))

                # 检查大小写问题
                if(input$kegg_species == "hsa") {
                  lower_case <- id_types$gene_symbols[grepl("^[a-z]", id_types$gene_symbols)]
                  if(length(lower_case) > 0) {
                    error_msg <- paste0(error_msg, "\n  大小写问题：", length(lower_case), "个基因是小写")
                    error_msg <- paste0(error_msg, "\n  建议：人类基因需要大写（如TP53，不是tp53）")
                  }
                } else if(input$kegg_species == "mmu") {
                  # 检查小鼠基因大小写
                  not_proper_case <- id_types$gene_symbols[!grepl("^[A-Z][a-z]+$", id_types$gene_symbols) & grepl("^[A-Za-z]", id_types$gene_symbols)]
                  if(length(not_proper_case) > 0) {
                    error_msg <- paste0(error_msg, "\n  大小写问题：", length(not_proper_case), "个基因大小写不正确")
                    error_msg <- paste0(error_msg, "\n  建议：小鼠基因需要首字母大写，其余小写（如Trp53，不是trp53或TRP53）")
                  }
                }
              }

              if(length(id_types$other_ids) > 0) {
                error_msg <- paste0(error_msg, "\n• 其他ID类型: ", length(id_types$other_ids), "个")
                error_msg <- paste0(error_msg, "\n  示例：", paste(head(id_types$other_ids, 3), collapse=", "))
                error_msg <- paste0(error_msg, "\n  建议：请检查这些ID的格式是否正确")
              }
            }

            showNotification(error_msg, type = "error", duration = 15)

            # 在调试模式下显示更多信息
            if(Sys.getenv("SHINY_DEBUG") == "TRUE" && !is.null(conversion_result$debug_info)) {
              cat("\n=== 背景基因转换详细调试信息 ===\n")
              print(conversion_result$debug_info)
            }
          }
        } else {
          showNotification(paste("错误：数据库包", db_pkg, "未安装"), type = "error")
        }
      }

      # ✅ 强制使用biofree.qyKEGGtools进行KEGG富集分析
      kegg_obj <- NULL

      if(!require("biofree.qyKEGGtools", quietly = TRUE)) {
        showNotification("❌ biofree.qyKEGGtools包未安装，无法进行KEGG富集分析", type = "error")
        return(NULL)
      }

      cat("✅ biofree.qyKEGGtools包已加载\n")
      cat(sprintf("📊 输入基因数量: %d\n", length(ids)))
      cat(sprintf("📊 输入基因ID类型: %s\n", class(ids)[1]))
      cat(sprintf("📊 输入基因示例: %s\n", paste(head(ids, 5), collapse = ", ")))

      # 🔧 转换为character向量
      ids_char <- as.character(ids)
      cat(sprintf("📊 转换为字符后示例: %s\n", paste(head(ids_char, 5), collapse = ", ")))

      # 构建参数列表（v2.1.0+ 支持 universe、pAdjustMethod、minGSSize、maxGSSize）
      args_list <- list(
        gene = ids_char,
        species = input$kegg_species,
        pCutoff = input$kegg_p,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500
      )
      if (!is.null(universe)) {
        universe_char <- as.character(universe)
        args_list$universe <- universe_char
        cat(sprintf("📊 背景基因数量: %d（已传入 enrich_local_KEGG）\n", length(universe_char)))
      }

      cat("🔧 开始调用biofree.qyKEGGtools::enrich_local_KEGG...\n")

      # 🔍 调试：检查函数是否存在
      cat(sprintf("📊 enrich_local_KEGG函数存在: %s\n", exists("enrich_local_KEGG", mode = "function")))

      # 🔍 调试：查看函数参数
      tryCatch({
        func_args <- formals(biofree.qyKEGGtools::enrich_local_KEGG)
        cat(sprintf("📊 enrich_local_KEGG函数参数: %s\n", paste(names(func_args), collapse = ", ")))
      }, error = function(e) {
        cat(sprintf("⚠️ 无法读取函数参数: %s\n", e$message))
      })

      # 使用do.call动态调用函数
      kegg_obj <- tryCatch({
        result <- do.call(biofree.qyKEGGtools::enrich_local_KEGG, args_list)
        cat("✅ biofree.qyKEGGtools::enrich_local_KEGG调用成功\n")

        # 🔍 检查返回结果
        if(!is.null(result)) {
          cat(sprintf("📊 结果对象类型: %s\n", class(result)[1]))
          if(inherits(result, "enrichResult")) {
            cat(sprintf("📊 富集结果数量: %d\n", nrow(result@result)))
          }
        }

        result
      }, error = function(e) {
        error_msg <- sprintf("biofree.qyKEGGtools调用失败: %s", e$message)
        cat(sprintf("⚠️ %s\n", error_msg))
        showNotification(error_msg, type = "error")
        NULL
      })

      # 🔧 如果返回NULL，给出详细错误信息
      if(is.null(kegg_obj)) {
        showNotification("❌ KEGG富集分析失败：biofree.qyKEGGtools返回NULL\n可能原因：基因ID格式不正确或本地KEGG数据库缺少对应物种的注释", type = "error")
        return(NULL)
      }

      # 检查结果
      result_df <- if(inherits(kegg_obj, "enrichResult")) {
        kegg_obj@result
      } else if(is.data.frame(kegg_obj)) {
        kegg_obj
      } else {
        cat(sprintf("⚠️ 未知的结果类型: %s\n", class(kegg_obj)[1]))
        showNotification("❌ KEGG富集分析返回了未知格式的结果", type = "error")
        return(NULL)
      }

      if(nrow(result_df) == 0) {
        showNotification("⚠️ KEGG富集分析没有找到显著富集的通路", type = "warning")
        return(NULL)
      }

      df <- result_df

      df$Description <- gsub(" - Mus musculus.*| - Homo sapiens.*", "", df$Description)

      db_pkg <- if(input$kegg_species == "mmu") "org.Mm.eg.db" else "org.Hs.eg.db"
      if(require(db_pkg, character.only = TRUE)) {
        db_obj <- get(db_pkg)
        all_entrez <- unique(unlist(strsplit(df$geneID, "/")))
        mapped <- AnnotationDbi::mapIds(db_obj, keys = all_entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

        df$geneID <- sapply(df$geneID, function(x) {
          ids <- unlist(strsplit(x, "/"))
          syms <- mapped[ids]
          syms[is.na(syms)] <- ids[is.na(syms)]
          paste(syms, collapse = "/")
        })
      }

      return(df)
    }, error = function(e) { showNotification(paste("KEGG Error:", e$message), type="error"); return(NULL) })
  })

  output$download_kegg <- downloadHandler(
    filename = function() {
      paste0("KEGG_Enrichment_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(kegg_data_processed())
      write.csv(kegg_data_processed(), file, row.names = FALSE)
    }
  )

  # KEGG图表生成reactive
  kegg_plot_obj <- reactive({
    req(kegg_data_processed())
    df <- kegg_data_processed()

    # 计算Fold Enrichment（如果数据框中没有）
    if (!"FoldEnrichment" %in% colnames(df)) {
      # 从GeneRatio和BgRatio计算Fold Enrichment
      # GeneRatio格式: "5/120", BgRatio格式: "50/5000"
      df$FoldEnrichment <- sapply(1:nrow(df), function(i) {
        gene_ratio <- as.numeric(strsplit(df$GeneRatio[i], "/")[[1]])
        bg_ratio <- as.numeric(strsplit(df$BgRatio[i], "/")[[1]])

        # 处理可能的NA或Inf值
        if (length(gene_ratio) < 2 || length(bg_ratio) < 2) {
          return(NA)
        }
        if (gene_ratio[2] == 0 || bg_ratio[2] == 0) {
          return(NA)
        }

        fe <- (gene_ratio[1] / gene_ratio[2]) / (bg_ratio[1] / bg_ratio[2])
        return(ifelse(is.finite(fe), fe, NA))
      })

      # 确保是数值类型
      df$FoldEnrichment <- as.numeric(df$FoldEnrichment)
    }

    df_plot <- head(df[order(df$p.adjust),], 20)

    txt_col <- if(input$theme_toggle) "white" else "black"
    grid_col <- if(input$theme_toggle) "#444444" else "#cccccc"

    font_face <- if(input$kegg_bold) "bold" else "plain"

    # 根据用户选择设置X轴变量
    if (input$kegg_x_axis == "FoldEnrichment") {
      x_var <- df_plot$FoldEnrichment
      x_label <- "Fold Enrichment"
    } else {
      x_var <- df_plot$Count
      x_label <- "Gene Count"
    }

    # 使用aes()而不是aes_string()
    p <- ggplot(df_plot, aes(x = x_var, y = reorder(Description, x_var), size = x_var, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = input$kegg_high_col, high = input$kegg_low_col) +
      theme_minimal() +
      labs(x = x_label, y = "", title = paste("KEGG Enrichment (", input$kegg_direction, ")")) +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        text = element_text(color = txt_col, size = input$kegg_font_size, face = font_face),
        axis.text = element_text(color = txt_col, size = input$kegg_font_size),
        legend.text = element_text(color = txt_col),
        legend.title = element_text(color = txt_col),
        axis.line = element_line(color = txt_col),
        panel.grid.major = element_line(color = grid_col),
        panel.grid.minor = element_line(color = grid_col)
      )

    return(p)
  })

  # KEGG图表下载
  output$download_kegg_plot <- downloadHandler(
    filename = function() {
      paste0("KEGG_Dotplot_", Sys.Date(), ".", input$kegg_export_format)
    },
    content = function(file) {
      req(kegg_plot_obj())

      # 获取当前图表
      p <- kegg_plot_obj()

      # 根据格式保存
      if (input$kegg_export_format == "png") {
        png(file, width = 10, height = 8, units = "in", res = 300)
      } else if (input$kegg_export_format == "pdf") {
        pdf(file, width = 10, height = 8)
      } else if (input$kegg_export_format == "svg") {
        svg(file, width = 10, height = 8)
      }

      print(p)
      dev.off()
    }
  )

  output$kegg_dotplot <- renderPlot({
    kegg_plot_obj()
  })

  output$kegg_table <- DT::renderDataTable({
    req(kegg_data_processed())
    DT::datatable(kegg_data_processed(), options = list(scrollX=T), rownames=F)
  })

  # =====================================================
  # 新增：单列基因 KEGG 富集分析
  # =====================================================

  # --- 单列基因 KEGG 富集分析 ---
  single_gene_kegg_data <- eventReactive(input$run_single_gene_kegg, {
    req(input$single_gene_file)

    showNotification("正在处理单列基因 KEGG 富集分析...", type = "message")

    # 读取单列基因文件
    gene_df <- read.csv(input$single_gene_file$datapath, header = TRUE)

    # 检查文件格式
    if (!"SYMBOL" %in% colnames(gene_df)) {
      showNotification("错误：CSV 文件必须包含 'SYMBOL' 列", type = "error")
      return(NULL)
    }

    # 获取基因列表
    gene_symbols <- gene_df$SYMBOL
    gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

    if (length(gene_symbols) == 0) {
      showNotification("错误：未找到有效的基因符号", type = "error")
      return(NULL)
    }

    showNotification(paste("找到", length(gene_symbols), "个基因进行富集分析"), type = "message")

    # 转换为 ENTREZID
    db_pkg <- if(input$single_gene_species == "mmu") "org.Mm.eg.db" else "org.Hs.eg.db"

    if (!require(db_pkg, character.only = TRUE, quietly = TRUE)) {
      showNotification(paste("错误：未安装", db_pkg), type = "error")
      return(NULL)
    }

    db_obj <- get(db_pkg)

    tryCatch({
      # 清理基因符号
      cleaned_genes <- clean_gene_symbols(gene_symbols, input$single_gene_species)

      # 使用智能转换函数将基因符号转换为ENTREZID
      conversion_result <- smart_gene_conversion(cleaned_genes, db_obj, "ENTREZID")

      if(is.null(conversion_result)) {
        # 如果智能转换失败，尝试直接使用SYMBOL keytype
        tryCatch({
          mapped <- AnnotationDbi::mapIds(db_obj,
                                         keys = cleaned_genes,
                                         column = "ENTREZID",
                                         keytype = "SYMBOL",
                                         multiVals = "first")
          entrez_ids <- na.omit(mapped)
        }, error = function(e) {
          showNotification(paste("基因符号转换失败:", e$message), type = "error")
          return(NULL)
        })
      } else {
        entrez_ids <- conversion_result$converted
        showNotification(paste("成功转换", length(entrez_ids), "个基因ID（通过",
                              conversion_result$keytype_used, "转换）"), type = "message")
      }

      if (length(entrez_ids) == 0) {
        showNotification("错误：无法将基因符号转换为 ENTREZID", type = "error")
        return(NULL)
      }

      # 对于单列基因分析，无法获取实验检测的背景基因集
      # 提供选项让用户选择是否使用全基因组作为背景
      universe <- NULL
      if(input$single_gene_kegg_use_background == "yes") {
        # 使用全基因组作为背景（虽然不是理想的选择，但比没有好）
        all_genes <- keys(db_obj, keytype = "ENTREZID")
        universe <- all_genes
        showNotification(paste("使用全基因组（", length(universe), "个基因）作为KEGG分析背景基因集"), type = "message")
      } else {
        showNotification("警告：未使用背景基因集，KEGG结果可能包含假阳性", type = "warning")
      }

      # 运行 KEGG 富集分析
      # 使用与主KEGG分析相同的逻辑，支持背景基因集
      kegg_obj <- NULL

      # 首先尝试使用 biofree.qyKEGGtools（v2.1.0+ 支持 universe）
      if(require("biofree.qyKEGGtools", quietly = TRUE)) {
        args_list <- list(
          gene = entrez_ids,
          species = input$single_gene_species,
          pCutoff = input$single_gene_kegg_p,
          pAdjustMethod = "BH",
          minGSSize = 10,
          maxGSSize = 500
        )
        if (!is.null(universe)) {
          args_list$universe <- as.character(universe)
        }
        kegg_obj <- tryCatch({
          do.call(biofree.qyKEGGtools::enrich_local_KEGG, args_list)
        }, error = function(e) {
          showNotification(paste("KEGG 分析出错:", e$message), type = "error")
          NULL
        })
      }

      # 如果biofree.qyKEGGtools不可用或失败，尝试使用clusterProfiler::enrichKEGG
      if(is.null(kegg_obj) && require("clusterProfiler", quietly = TRUE)) {
        showNotification("使用clusterProfiler::enrichKEGG进行单列基因KEGG富集分析", type = "message")

        # 设置KEGG数据库的物种代码
        kegg_org <- if(input$single_gene_species == "mmu") "mmu" else "hsa"

        kegg_obj <- clusterProfiler::enrichKEGG(
          gene = entrez_ids,
          organism = kegg_org,
          pvalueCutoff = input$single_gene_kegg_p,
          pAdjustMethod = "BH",
          universe = universe
        )
      }

      if (is.null(kegg_obj) || nrow(kegg_obj@result) == 0) {
        showNotification("KEGG富集分析没有结果。", type = "warning")
        return(NULL)
      }

      df <- kegg_obj@result

      # 清理描述信息
      df$Description <- gsub(" - Mus musculus.*| - Homo sapiens.*", "", df$Description)

      # 将 ENTREZID 转换回 SYMBOL 用于显示
      all_entrez <- unique(unlist(strsplit(df$geneID, "/")))
      symbol_mapped <- AnnotationDbi::mapIds(db_obj,
                                             keys = all_entrez,
                                             column = "SYMBOL",
                                             keytype = "ENTREZID",
                                             multiVals = "first")

      df$geneID <- sapply(df$geneID, function(x) {
        ids <- unlist(strsplit(x, "/"))
        syms <- symbol_mapped[ids]
        syms[is.na(syms)] <- ids[is.na(syms)]
        paste(syms, collapse = "/")
      })

      # 添加原始基因数量信息
      attr(df, "input_genes_count") <- length(gene_symbols)
      attr(df, "mapped_genes_count") <- length(entrez_ids)

      return(df)

    }, error = function(e) {
      showNotification(paste("KEGG 分析错误:", e$message), type = "error")
      return(NULL)
    })
  })

  # --- 单列基因 KEGG 结果下载 ---
  output$download_single_gene_kegg <- downloadHandler(
    filename = function() {
      paste0("Single_Gene_KEGG_Enrichment_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(single_gene_kegg_data())
      write.csv(single_gene_kegg_data(), file, row.names = FALSE)
    }
  )

  # 单列基因KEGG图表生成reactive
  single_gene_kegg_plot_obj <- reactive({
    req(single_gene_kegg_data())
    df <- single_gene_kegg_data()

    # 计算Fold Enrichment（如果数据框中没有）
    if (!"FoldEnrichment" %in% colnames(df)) {
      # 从GeneRatio和BgRatio计算Fold Enrichment
      # GeneRatio格式: "5/120", BgRatio格式: "50/5000"
      df$FoldEnrichment <- sapply(1:nrow(df), function(i) {
        gene_ratio <- as.numeric(strsplit(df$GeneRatio[i], "/")[[1]])
        bg_ratio <- as.numeric(strsplit(df$BgRatio[i], "/")[[1]])

        # 处理可能的NA或Inf值
        if (length(gene_ratio) < 2 || length(bg_ratio) < 2) {
          return(NA)
        }
        if (gene_ratio[2] == 0 || bg_ratio[2] == 0) {
          return(NA)
        }

        fe <- (gene_ratio[1] / gene_ratio[2]) / (bg_ratio[1] / bg_ratio[2])
        return(ifelse(is.finite(fe), fe, NA))
      })

      # 确保是数值类型
      df$FoldEnrichment <- as.numeric(df$FoldEnrichment)
    }

    df_plot <- head(df[order(df$p.adjust),], 20)

    input_count <- attr(df, "input_genes_count")
    mapped_count <- attr(df, "mapped_genes_count")

    txt_col <- if(input$theme_toggle) "white" else "black"
    grid_col <- if(input$theme_toggle) "#444444" else "#cccccc"

    font_face <- if(input$single_gene_kegg_bold) "bold" else "plain"

    # 根据用户选择设置X轴变量
    if (input$single_gene_kegg_x_axis == "FoldEnrichment") {
      x_var <- df_plot$FoldEnrichment
      x_label <- "Fold Enrichment"
    } else {
      x_var <- df_plot$Count
      x_label <- "Gene Count"
    }

    # 使用aes()而不是aes_string()
    p <- ggplot(df_plot, aes(x = x_var, y = reorder(Description, x_var), size = x_var, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = input$single_gene_kegg_high_col, high = input$single_gene_kegg_low_col) +
      theme_minimal() +
      labs(
        x = x_label,
        y = "",
        title = "单列基因 KEGG 富集分析",
        subtitle = paste("输入基因:", input_count, "个 | 成功映射:", mapped_count, "个")
      ) +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = txt_col, hjust = 0.5),
        text = element_text(color = txt_col, size = input$single_gene_kegg_font_size, face = font_face),
        axis.text = element_text(color = txt_col, size = input$single_gene_kegg_font_size),
        legend.text = element_text(color = txt_col),
        legend.title = element_text(color = txt_col),
        axis.line = element_line(color = txt_col),
        panel.grid.major = element_line(color = grid_col),
        panel.grid.minor = element_line(color = grid_col)
      )

    return(p)
  })

  # 单列基因KEGG图表下载
  output$download_single_gene_kegg_plot <- downloadHandler(
    filename = function() {
      paste0("Single_Gene_KEGG_Dotplot_", Sys.Date(), ".", input$single_gene_kegg_export_format)
    },
    content = function(file) {
      req(single_gene_kegg_plot_obj())

      # 获取当前图表
      p <- single_gene_kegg_plot_obj()

      # 根据格式保存
      if (input$single_gene_kegg_export_format == "png") {
        png(file, width = 10, height = 8, units = "in", res = 300)
      } else if (input$single_gene_kegg_export_format == "pdf") {
        pdf(file, width = 10, height = 8)
      } else if (input$single_gene_kegg_export_format == "svg") {
        svg(file, width = 10, height = 8)
      }

      print(p)
      dev.off()
    }
  )

  # --- 单列基因 KEGG 点图 ---
  output$single_gene_kegg_dotplot <- renderPlot({
    single_gene_kegg_plot_obj()
  })

  # --- 单列基因 KEGG 结果表格 ---
  output$single_gene_kegg_table <- DT::renderDataTable({
    req(single_gene_kegg_data())
    DT::datatable(single_gene_kegg_data(), options = list(scrollX = TRUE), rownames = FALSE)
  })

  # 返回KEGG结果供其他模块使用
  return(reactive({ kegg_data_processed() }))
}