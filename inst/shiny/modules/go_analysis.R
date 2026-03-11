# =====================================================
# GO富集分析模块
# =====================================================

go_analysis_server <- function(input, output, session, deg_results) {

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
  # GO 分析 - 基于差异基因
  # =====================================================
  go_data_processed <- eventReactive(input$run_go, {
    req(deg_results())

    # 从deg_results中提取差异分析结果和背景基因
    deg_data <- deg_results()
    res_df <- deg_data$deg_df
    background_genes <- deg_data$background_genes

    # 清理背景基因符号
    if (!is.null(background_genes) && length(background_genes) > 0) {
      background_genes <- clean_gene_symbols(background_genes, input$go_species)
    }

    # 显示进度提示
    withProgress(message = 'GO分析进行中...', value = 0, {
      # 进度1: 准备数据 (10%)
      incProgress(0.1, detail = "准备分析数据...")

      target_status <- switch(input$go_direction, "Up" = "Up", "Down" = "Down", "All" = c("Up", "Down"))
      ids <- res_df %>% dplyr::filter(Status %in% target_status & !is.na(ENTREZID)) %>% dplyr::pull(ENTREZID)

      if(length(ids) == 0) {
        showNotification("无有效ENTREZID，请检查基因注释结果", type="error")
        return(NULL)
      }

      tryCatch({
        # 进度2: 加载数据库 (20%)
        incProgress(0.1, detail = "加载基因注释数据库...")

        # 根据物种选择对应的org包
        db_pkg <- if(input$go_species == "mmu") "org.Mm.eg.db" else "org.Hs.eg.db"
        if(!require(db_pkg, character.only = TRUE, quietly = TRUE)) {
          showNotification(paste("请先安装", db_pkg, "包"), type="error")
          return(NULL)
        }

        db_obj <- get(db_pkg)

        # 准备背景基因集（如果可用）
        universe <- NULL
        if(!is.null(background_genes) && length(background_genes) > 0) {
          # 清理背景基因符号
          cleaned_background <- clean_gene_symbols(background_genes, input$go_species)

          # 将背景基因符号转换为ENTREZID（使用智能转换）
          conversion_result <- smart_gene_conversion(cleaned_background, db_obj, "ENTREZID")

          # 检查转换结果
          if(!is.null(conversion_result$converted) && length(conversion_result$converted) > 0) {
            bg_entrez <- conversion_result$converted
            universe <- bg_entrez

            # 显示成功信息
            if(!is.null(conversion_result$keytype_used)) {
              showNotification(paste("使用", length(universe), "个检测到的基因作为背景基因集（通过",
                                    conversion_result$keytype_used, "转换）"), type = "message")
            } else {
              showNotification(paste("使用", length(universe), "个检测到的基因作为背景基因集"), type = "message")
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
              id_types <- identify_gene_id_types(sample_genes, input$go_species)
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
                if(input$go_species == "hsa") {
                  lower_case <- id_types$gene_symbols[grepl("^[a-z]", id_types$gene_symbols)]
                  if(length(lower_case) > 0) {
                    error_msg <- paste0(error_msg, "\n  大小写问题：", length(lower_case), "个基因是小写")
                    error_msg <- paste0(error_msg, "\n  建议：人类基因需要大写（如TP53，不是tp53）")
                  }
                } else if(input$go_species == "mmu") {
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
        }

        # 进度3: 进行GO富集分析 (40%)
        incProgress(0.2, detail = "进行GO富集分析...")

        # 进行GO富集分析（使用背景基因集）
        go_obj <- clusterProfiler::enrichGO(
          gene          = ids,
          OrgDb         = db_obj,
          keyType       = "ENTREZID",
          ont           = input$go_ontology,      # BP, MF, CC
          pAdjustMethod = "BH",
          pvalueCutoff  = input$go_p,
          qvalueCutoff  = 0.2,
          readable      = TRUE,                   # 将ENTREZID转换为Gene Symbol
          universe      = universe                # 关键改进：使用背景基因集
        )

        if(is.null(go_obj) || nrow(go_obj@result) == 0) {
          showNotification("GO富集分析没有结果。", type = "warning")
          return(NULL)
        }

        # 进度4: 处理结果数据 (60%)
        incProgress(0.2, detail = "处理分析结果...")

        df <- go_obj@result

        # 简化描述（移除物种信息）
        df$Description <- gsub("\\s*\\(.*\\)$", "", df$Description)

        # 将GeneRatio从字符格式转换为数值比例
        if ("GeneRatio" %in% colnames(df)) {
          df$GeneRatio_numeric <- sapply(strsplit(df$GeneRatio, "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
        }

        # 将BgRatio从字符格式转换为数值比例
        if ("BgRatio" %in% colnames(df)) {
          df$BgRatio_numeric <- sapply(strsplit(df$BgRatio, "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
        }

        # 计算富集倍数（Fold Enrichment）
        if ("GeneRatio_numeric" %in% colnames(df) && "BgRatio_numeric" %in% colnames(df)) {
          df$Fold_Enrichment <- df$GeneRatio_numeric / df$BgRatio_numeric

          # 处理除零错误（如果BgRatio为0）
          df$Fold_Enrichment[is.infinite(df$Fold_Enrichment) | is.na(df$Fold_Enrichment)] <- NA
        }

        # 进度5: 完成 (100%)
        incProgress(0.4, detail = "分析完成！")
        Sys.sleep(0.5)  # 短暂延迟，让用户看到完成消息

        return(df)

      }, error = function(e) {
        showNotification(paste("GO分析错误:", e$message), type="error")
        return(NULL)
      })
    })
  })

  # GO结果下载
  output$download_go <- downloadHandler(
    filename = function() {
      paste0("GO_Enrichment_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- go_data_processed()
      if(!is.null(df)) {
        write.csv(df, file, row.names = FALSE)
      }
    }
  )

  # GO图表生成reactive
  go_plot_obj <- reactive({
    df <- go_data_processed()
    if(is.null(df)) return(NULL)

    # 根据用户选择确定X轴变量
    x_var <- switch(input$go_x_axis,
                    "GeneRatio" = "GeneRatio_numeric",
                    "BgRatio" = "BgRatio_numeric",
                    "Count" = "Count",
                    "FoldEnrichment" = "Fold_Enrichment")

    # 根据X轴变量排序数据，并选择Top N
    # 降序排列（ Fold Enrichment最大的在上面）
    df_sorted <- df[order(df[[x_var]], decreasing = TRUE), ]
    top_n <- min(input$go_top_n, nrow(df_sorted))
    df_top <- df_sorted[1:top_n, ]

    x_label <- switch(input$go_x_axis,
                      "GeneRatio" = "Gene Ratio",
                      "BgRatio" = "Background Ratio",
                      "Count" = "Gene Count",
                      "FoldEnrichment" = "Fold Enrichment")

    # 创建点图 - Y轴已经按X轴变量排序，所以直接使用Description即可
    p <- ggplot(df_top, aes(x = .data[[x_var]], y = reorder(Description, .data[[x_var]]))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = input$go_low_col, high = input$go_high_col) +
      theme_bw(base_size = input$go_font_size) +
      labs(x = x_label, y = "GO Term",
           title = paste("GO Enrichment Analysis (", input$go_ontology, ")", sep = ""),
           size = "Gene Count", color = "Adj. P-value") +
      theme(axis.text.y = element_text(face = ifelse(input$go_bold, "bold", "plain")))

    return(p)
  })

  # GO图表下载
  output$download_go_plot <- downloadHandler(
    filename = function() {
      paste0("GO_Dotplot_", Sys.Date(), ".", input$go_export_format)
    },
    content = function(file) {
      p <- go_plot_obj()
      if(is.null(p)) return(NULL)

      # 根据格式保存
      if (input$go_export_format == "png") {
        png(file, width = 10, height = 8, units = "in", res = 300)
      } else if (input$go_export_format == "pdf") {
        pdf(file, width = 10, height = 8)
      } else if (input$go_export_format == "svg") {
        svg(file, width = 10, height = 8)
      }

      print(p)
      dev.off()
    }
  )

  # GO点图
  output$go_dotplot <- renderPlot({
    go_plot_obj()
  })

  # GO结果表格
  output$go_table <- DT::renderDataTable({
    df <- go_data_processed()
    if(is.null(df)) return(NULL)

    # 选择显示的列（添加Fold_Enrichment）
    display_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "Fold_Enrichment",
                      "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    display_cols <- display_cols[display_cols %in% colnames(df)]

    # 格式化富集倍数（保留2位小数）
    if ("Fold_Enrichment" %in% colnames(df)) {
      df$Fold_Enrichment <- round(df$Fold_Enrichment, 2)
    }

    DT::datatable(df[, display_cols],
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf')
                  ),
                  rownames = FALSE,
                  class = "display compact")
  })

  # =====================================================
  # 单列基因 GO 富集分析
  # =====================================================
  single_gene_go_data <- eventReactive(input$run_single_gene_go, {
    req(input$single_gene_go_file)

    tryCatch({
      # 读取基因列表文件
      gene_df <- read.csv(input$single_gene_go_file$datapath)
      if(!"SYMBOL" %in% colnames(gene_df)) {
        showNotification("CSV文件必须包含'SYMBOL'列", type="error")
        return(NULL)
      }

      gene_symbols <- gene_df$SYMBOL

      # 根据物种选择对应的org包
      db_pkg <- if(input$single_gene_species == "mmu") "org.Mm.eg.db" else "org.Hs.eg.db"
      if(!require(db_pkg, character.only = TRUE, quietly = TRUE)) {
        showNotification(paste("请先安装", db_pkg, "包"), type="error")
        return(NULL)
      }

      # 将Gene Symbol转换为ENTREZID（使用智能转换）
      db_obj <- get(db_pkg)

      # 清理基因符号
      cleaned_genes <- clean_gene_symbols(gene_symbols, input$single_gene_species)

      # 使用智能转换函数
      conversion_result <- smart_gene_conversion(cleaned_genes, db_obj, "ENTREZID")

      if(is.null(conversion_result)) {
        # 如果智能转换失败，尝试直接使用SYMBOL keytype
        tryCatch({
          entrez_ids <- AnnotationDbi::mapIds(db_obj,
                                             keys = cleaned_genes,
                                             column = "ENTREZID",
                                             keytype = "SYMBOL",
                                             multiVals = "first")
          entrez_ids <- entrez_ids[!is.na(entrez_ids)]
        }, error = function(e) {
          showNotification(paste("基因符号转换失败:", e$message), type = "error")
          return(NULL)
        })
      } else {
        entrez_ids <- conversion_result$converted
        showNotification(paste("成功转换", length(entrez_ids), "个基因ID（通过",
                              conversion_result$keytype_used, "转换）"), type = "message")
      }

      if(length(entrez_ids) == 0) {
        showNotification("无法将基因符号转换为ENTREZID", type="error")
        return(NULL)
      }

      # 对于单列基因分析，无法获取背景基因集
      # 显示警告信息，说明这可能导致假阳性
      if(input$single_gene_use_background == "no") {
        universe <- NULL
        showNotification("警告：未使用背景基因集，结果可能包含假阳性", type = "warning")
      } else {
        # 尝试获取所有可能的基因作为背景（全基因组）
        # 注意：这仍然不是理想的背景，但比没有好
        all_genes <- keys(db_obj, keytype = "ENTREZID")
        universe <- all_genes
        showNotification(paste("使用全基因组（", length(universe), "个基因）作为背景基因集"), type = "message")
      }

      # 进行GO富集分析
      go_obj <- clusterProfiler::enrichGO(
        gene          = entrez_ids,
        OrgDb         = db_obj,
        keyType       = "ENTREZID",
        ont           = input$single_gene_go_ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = input$single_gene_go_p,
        qvalueCutoff  = 0.2,
        readable      = TRUE,
        universe      = universe  # 使用背景基因集
      )

      if(is.null(go_obj) || nrow(go_obj@result) == 0) {
        showNotification("单列基因GO富集分析没有结果。", type = "warning")
        return(NULL)
      }

      df <- go_obj@result
      df$Description <- gsub("\\s*\\(.*\\)$", "", df$Description)

      # 将GeneRatio从字符格式转换为数值比例
      if ("GeneRatio" %in% colnames(df)) {
        df$GeneRatio_numeric <- sapply(strsplit(df$GeneRatio, "/"), function(x) {
          as.numeric(x[1]) / as.numeric(x[2])
        })
      }

      # 将BgRatio从字符格式转换为数值比例
      if ("BgRatio" %in% colnames(df)) {
        df$BgRatio_numeric <- sapply(strsplit(df$BgRatio, "/"), function(x) {
          as.numeric(x[1]) / as.numeric(x[2])
        })
      }

      # 计算富集倍数（Fold Enrichment）
      if ("GeneRatio_numeric" %in% colnames(df) && "BgRatio_numeric" %in% colnames(df)) {
        df$Fold_Enrichment <- df$GeneRatio_numeric / df$BgRatio_numeric

        # 处理除零错误（如果BgRatio为0）
        df$Fold_Enrichment[is.infinite(df$Fold_Enrichment) | is.na(df$Fold_Enrichment)] <- NA
      }

      return(df)

    }, error = function(e) {
      showNotification(paste("单列基因GO分析错误:", e$message), type="error")
      return(NULL)
    })
  })

  # 单列基因GO结果下载
  output$download_single_gene_go <- downloadHandler(
    filename = function() {
      paste0("Single_Gene_GO_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- single_gene_go_data()
      if(!is.null(df)) {
        write.csv(df, file, row.names = FALSE)
      }
    }
  )

  # 单列基因GO图表生成reactive
  single_gene_go_plot_obj <- reactive({
    df <- single_gene_go_data()
    if(is.null(df)) return(NULL)

    # 根据用户选择确定X轴变量
    x_var <- switch(input$single_gene_go_x_axis,
                    "GeneRatio" = "GeneRatio_numeric",
                    "BgRatio" = "BgRatio_numeric",
                    "Count" = "Count",
                    "FoldEnrichment" = "Fold_Enrichment")

    # 根据X轴变量排序数据，并选择Top N
    # 降序排列（ Fold Enrichment最大的在上面）
    df_sorted <- df[order(df[[x_var]], decreasing = TRUE), ]
    top_n <- min(20, nrow(df_sorted))
    df_top <- df_sorted[1:top_n, ]

    x_label <- switch(input$single_gene_go_x_axis,
                      "GeneRatio" = "Gene Ratio",
                      "BgRatio" = "Background Ratio",
                      "Count" = "Gene Count",
                      "FoldEnrichment" = "Fold Enrichment")

    # 创建点图
    p <- ggplot(df_top, aes(x = .data[[x_var]], y = reorder(Description, .data[[x_var]]))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = input$single_gene_go_low_col, high = input$single_gene_go_high_col) +
      theme_bw(base_size = input$single_gene_go_font_size) +
      labs(x = x_label, y = "GO Term",
           title = paste("Single Gene GO Enrichment (", input$single_gene_go_ontology, ")", sep = ""),
           size = "Gene Count", color = "Adj. P-value") +
      theme(axis.text.y = element_text(face = ifelse(input$single_gene_go_bold, "bold", "plain")))

    return(p)
  })

  # 单列基因GO图表下载
  output$download_single_gene_go_plot <- downloadHandler(
    filename = function() {
      paste0("Single_Gene_GO_Dotplot_", Sys.Date(), ".", input$single_gene_go_export_format)
    },
    content = function(file) {
      p <- single_gene_go_plot_obj()
      if(is.null(p)) return(NULL)

      # 根据格式保存
      if (input$single_gene_go_export_format == "png") {
        png(file, width = 10, height = 8, units = "in", res = 300)
      } else if (input$single_gene_go_export_format == "pdf") {
        pdf(file, width = 10, height = 8)
      } else if (input$single_gene_go_export_format == "svg") {
        svg(file, width = 10, height = 8)
      }

      print(p)
      dev.off()
    }
  )

  # 单列基因GO点图
  output$single_gene_go_dotplot <- renderPlot({
    single_gene_go_plot_obj()
  })

  # 单列基因GO结果表格
  output$single_gene_go_table <- DT::renderDataTable({
    df <- single_gene_go_data()
    if(is.null(df)) return(NULL)

    # 选择显示的列（添加Fold_Enrichment）
    display_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "Fold_Enrichment",
                      "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    display_cols <- display_cols[display_cols %in% colnames(df)]

    # 格式化富集倍数（保留2位小数）
    if ("Fold_Enrichment" %in% colnames(df)) {
      df$Fold_Enrichment <- round(df$Fold_Enrichment, 2)
    }

    DT::datatable(df[, display_cols],
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf')
                  ),
                  rownames = FALSE,
                  class = "display compact")
  })

  # 返回GO结果供其他模块使用
  return(reactive({ go_data_processed() }))
}