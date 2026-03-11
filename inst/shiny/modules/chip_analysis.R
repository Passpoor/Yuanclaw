# =====================================================
# 芯片数据分析模块 (Microarray Data Analysis)
# =====================================================
# 功能：
# 1. 解析 GEO Series Matrix 文件
# 2. 解析 SOFT 平台注释文件
# 3. 探针注释 (探针ID → 基因符号)
# 4. limma 差异分析
# 5. 与现有下游模块集成 (KEGG/GO/GSEA/TF/通路活性)
# =====================================================

# =====================================================
# 1. GEO Series Matrix 文件解析
# =====================================================

#' 解析 GEO Series Matrix 文件
#'
#' @param file_path Series Matrix 文件路径
#' @return list 包含表达矩阵和元数据
parse_geo_series_matrix <- function(file_path) {
  cat("📂 开始解析 GEO Series Matrix 文件...\n")

  tryCatch({
    # 读取所有行
    lines <- readLines(file_path, warn = FALSE)

    # 查找数据起始标记
    start_idx <- which(lines == "!series_matrix_table_begin")

    if (length(start_idx) == 0) {
      return(list(
        success = FALSE,
        error = "未找到 !series_matrix_table_begin 标记，这不是有效的 GEO Series Matrix 文件"
      ))
    }

    start_idx <- start_idx[1] + 1  # 跳过标记行本身

    # 查找数据结束标记（可选）
    end_idx <- which(lines == "!series_matrix_table_end")
    if (length(end_idx) > 0) {
      end_idx <- end_idx[1] - 1
    } else {
      # 如果没有结束标记，读取到文件末尾
      end_idx <- length(lines)
    }

    cat(sprintf("✅ 找到数据区域: 第 %d - %d 行\n", start_idx, end_idx))

    # 提取矩阵数据
    matrix_lines <- lines[start_idx:end_idx]

    # 读取为数据框
    text_connection <- textConnection(matrix_lines)
    expr_matrix <- read.table(
      text_connection,
      header = TRUE,
      row.names = 1,
      sep = "\t",
      quote = "\"",  # 允许引号
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    close(text_connection)

    # 去除列名中的引号
    colnames(expr_matrix) <- gsub('"', '', colnames(expr_matrix))
    rownames(expr_matrix) <- gsub('"', '', rownames(expr_matrix))

    cat(sprintf("✅ 清理引号后样本名示例: %s\n", paste(head(colnames(expr_matrix), 3), collapse = ", ")))
    cat(sprintf("✅ 清理引号后探针ID示例: %s\n", paste(head(rownames(expr_matrix), 3), collapse = ", ")))

    # 提取元数据（用于样本分组提示）
    metadata <- extract_geo_metadata(lines[1:(start_idx-2)])

    cat(sprintf("✅ 解析完成: %d 探针 × %d 样本\n",
                nrow(expr_matrix), ncol(expr_matrix)))

    return(list(
      success = TRUE,
      matrix = expr_matrix,
      metadata = metadata,
      n_probes = nrow(expr_matrix),
      n_samples = ncol(expr_matrix),
      sample_names = colnames(expr_matrix)
    ))

  }, error = function(e) {
    return(list(
      success = FALSE,
      error = paste("解析文件时出错:", e$message)
    ))
  })
}

#' 提取 GEO 元数据
#'
#' @param metadata_lines 元数据行
#' @return list 包含样本描述等信息
extract_geo_metadata <- function(metadata_lines) {
  result <- list(
    has_metadata = FALSE,
    sample_descriptions = NULL,
    sample_titles = NULL
  )

  # 查找样本描述行
  desc_line <- grep("^!Sample_description", metadata_lines, value = TRUE)

  if (length(desc_line) > 0) {
    # 提取描述信息
    descriptions <- gsub("^!Sample_description\t", "", desc_line)
    descriptions <- gsub('"', '', descriptions)
    # 分割成向量
    descriptions <- strsplit(descriptions[1], "\t")[[1]]
    result$sample_descriptions <- descriptions
    result$has_metadata <- TRUE
  }

  # 查找样本标题行
  title_line <- grep("^!Sample_title", metadata_lines, value = TRUE)

  if (length(title_line) > 0) {
    # 提取标题信息
    titles <- gsub("^!Sample_title\t", "", title_line)
    titles <- gsub('"', '', titles)
    # 分割成向量
    titles <- strsplit(titles[1], "\t")[[1]]
    result$sample_titles <- titles
    result$has_metadata <- TRUE
  }

  return(result)
}

# =====================================================
# 2. SOFT 平台文件解析
# =====================================================

#' 解析 SOFT 平台注释文件
#'
#' @param file_path SOFT 文件路径
#' @param separator 字段分隔符（正则表达式，默认为Tab）
#' @return list 包含原始表格和探针-基因映射
parse_platform_annotation <- function(file_path, separator = "\t") {
  cat("📋 开始解析 SOFT 平台文件...\n")
  cat(sprintf("📋 使用分隔符: %s\n", separator))

  tryCatch({
    lines <- readLines(file_path, warn = FALSE)

    # 查找平台表起始标记
    start_idx <- which(lines == "!platform_table_begin")

    if (length(start_idx) == 0) {
      return(list(
        success = FALSE,
        error = "未找到 !platform_table_begin 标记"
      ))
    }

    start_idx <- start_idx[1] + 1

    # 读取到文件末尾或下一个 ^ 标记
    remaining_lines <- lines[start_idx:length(lines)]
    next_section <- which(grepl("^\\^", remaining_lines))

    if (length(next_section) > 0) {
      end_idx <- start_idx + next_section[1] - 2
    } else {
      end_idx <- length(lines)
    }

    table_lines <- lines[start_idx:end_idx]

    # 读取表头（使用用户指定的分隔符）
    header <- strsplit(table_lines[1], separator)[[1]]

    # 读取数据（使用用户指定的分隔符）
    text_conn <- textConnection(table_lines)
    raw_table <- read.table(
      text_conn,
      header = TRUE,
      sep = separator,
      quote = "",
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE
    )
    close(text_conn)

    cat(sprintf("📋 平台注释文件解析: %d 列 × %d 行\n",
                ncol(raw_table), nrow(raw_table)))
    cat("📋 列名:", paste(colnames(raw_table), collapse = ", "), "\n")

    # 🔍 显示各列的示例数据，帮助用户选择
    cat("📋 各列示例数据:\n")
    for (col in colnames(raw_table)) {
      sample_vals <- head(raw_table[[col]][!is.na(raw_table[[col]]) & raw_table[[col]] != ""], 3)
      cat(sprintf("   %s: %s\n", col, paste(sample_vals, collapse = ", ")))
    }

    # ❌ 不再自动检测基因符号列，由用户手动选择
    cat("⚠️ 请用户手动选择基因符号列\n")

    # 返回原始表格，不进行自动检测和映射
    return(list(
      success = TRUE,
      raw_table = raw_table,
      mapping = NULL,  # 不提供自动映射
      gene_symbol_col = NULL,  # 不提供自动检测的列名
      needs_manual_selection = TRUE,  # 标记需要手动选择
      message = "请手动选择ID列和基因列"
    ))

  }, error = function(e) {
    return(list(
      success = FALSE,
      error = paste("解析 SOFT 文件时出错:", e$message)
    ))
  })
}

#' 智能检测基因符号列（改进版 - 基于内容分析）
#'
#' @param table 平台注释数据框
#' @return 字符列名或 NULL
detect_gene_symbol_column <- function(table) {
  cat("🔍 开始智能检测基因符号列...\n")

  # 常见的基因符号列名（优先级高）
  high_priority_names <- c(
    "GENE_SYMBOL",
    "Gene.Symbol",
    "Gene_Symbol",
    "gene_symbol",
    "SYMBOL",
    "Symbol",
    "Gene Symbol"
  )

  # 方法1: 优先级列名匹配
  for (name in high_priority_names) {
    if (name %in% colnames(table)) {
      # 验证：检查该列是否真的包含基因符号格式
      col_data <- table[[name]]
      col_data <- col_data[!is.na(col_data) & col_data != ""]
      n_check <- min(50, length(col_data))

      if (n_check >= 10) {
        # 检查是否包含典型基因符号格式
        # 基因符号：字母开头，可能包含数字和连字符，长度2-20
        pattern <- "^[A-Z][A-Z0-9\\-]{1,15}$"
        match_ratio <- sum(grepl(pattern, col_data[1:n_check])) / n_check

        if (match_ratio > 0.5) {  # 超过50%匹配
          cat(sprintf("✅ 智能检测到基因列(高优先级): %s (匹配率: %.1f%%)\n",
                      name, match_ratio*100))
          return(name)
        }
      }
    }
  }

  # 方法2: 列名部分匹配
  for (name in high_priority_names) {
    matching_cols <- colnames(table)[sapply(colnames(table), function(col) {
      grepl(name, col, ignore.case = TRUE)
    })]

    if (length(matching_cols) > 0) {
      # 验证内容
      for (col_name in matching_cols) {
        col_data <- table[[col_name]]
        col_data <- col_data[!is.na(col_data) & col_data != ""]
        n_check <- min(50, length(col_data))

        if (n_check >= 10) {
          pattern <- "^[A-Z][A-Z0-9\\-]{1,15}$"
          match_ratio <- sum(grepl(pattern, col_data[1:n_check])) / n_check

          if (match_ratio > 0.5) {
            cat(sprintf("✅ 智能检测到基因列(模糊匹配): %s (匹配率: %.1f%%)\n",
                        col_name, match_ratio*100))
            return(col_name)
          }
        }
      }
    }
  }

  # 方法3: 检查所有列的内容特征
  best_col <- NULL
  best_match_ratio <- 0

  for (col_name in colnames(table)) {
    # 跳过明显不是基因的列
    if (col_name %in% c("ID", "SPOT_ID", "CONTROL_TYPE", "REFSEQ", "GB_ACC",
                       "UNIGENE_ID", "ENSEMBL_ID", "TIGR_ID", "ACCESSION_STRING",
                       "CHROMOSOMAL_LOCATION", "CYTOBAND", "DESCRIPTION", "GO_ID",
                       "SEQUENCE")) {
      next
    }

    col_data <- table[[col_name]]
    col_data <- col_data[!is.na(col_data) & col_data != ""]
    n_check <- min(100, length(col_data))

    if (n_check < 10) next

    # 检查基因符号特征
    # 典型基因符号：TP53, EGFR, BRCA1, MYC, IL-6, TNF-alpha
    patterns <- c(
      "^[A-Z][A-Z0-9]{1,10}$",      # 简单基因符号：TP53, EGFR
      "^[A-Z][A-Z0-9]{1,5}-[0-9]+$", # 带数字的：BRCA1-001
      "^[A-Z]{2,6}-[A-Z0-9]{1,3}$"   # 带连字符的：IL-6, TNF-a
    )

    match_count <- sum(sapply(patterns, function(p) {
      sum(grepl(p, col_data[1:n_check]))
    }))

    match_ratio <- match_count / n_check

    if (match_ratio > best_match_ratio && match_ratio > 0.3) {
      best_match_ratio <- match_ratio
      best_col <- col_name
    }
  }

  if (!is.null(best_col)) {
    cat(sprintf("✅ 智能检测到基因列(内容分析): %s (匹配率: %.1f%%)\n",
                best_col, best_match_ratio*100))
    return(best_col)
  }

  cat("⚠️ 无法自动检测基因符号列，需要用户手动选择\n")
  return(NULL)
}

#' 执行探针注释（合并表达矩阵和基因映射）
#'
#' @param expr_matrix 探针表达矩阵
#' @param probe_gene_map 探针-基因映射数据框
#' @return 基因表达矩阵
annotate_probe_matrix <- function(expr_matrix, probe_gene_map) {
  cat("🔄 开始探针注释...\n")

  # 找到共同的探针
  common_probes <- intersect(rownames(expr_matrix), probe_gene_map$probe_id)

  if (length(common_probes) == 0) {
    stop("❌ 探针ID不匹配，请检查平台注释文件是否正确")
  }

  cat(sprintf("📊 共同探针数: %d (表达矩阵: %d, 注释文件: %d)\n",
              length(common_probes),
              nrow(expr_matrix),
              nrow(probe_gene_map)))

  # 提取共同探针的表达数据
  expr_subset <- expr_matrix[common_probes, , drop = FALSE]

  # 添加探针ID列以便合并
  expr_subset <- data.frame(probe_id = rownames(expr_subset), expr_subset,
                           stringsAsFactors = FALSE)

  # 合并基因符号
  merged <- merge(expr_subset, probe_gene_map, by = "probe_id")

  # 移除探针ID列
  probe_id_col <- which(names(merged) == "probe_id")
  if (length(probe_id_col) > 0) {
    merged <- merged[, -probe_id_col]
  }

  # 按基因分组，对表达值取平均值
  # 如果一个基因对应多个探针，取平均值
  library(dplyr)

  expr_annotated <- merged %>%
    group_by(gene_symbol) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames("gene_symbol")

  # 转换为矩阵
  expr_annotated <- as.matrix(expr_annotated)

  cat(sprintf("✅ 探针注释完成: %d 探针 → %d 基因\n",
              length(common_probes),
              nrow(expr_annotated)))

  # 检查是否有基因被重复注释
  n_probes_per_gene <- table(merged$gene_symbol)
  n_multi_probe_genes <- sum(n_probes_per_gene > 1)

  if (n_multi_probe_genes > 0) {
    cat(sprintf("📊 其中 %d 个基因有多个探针(已取平均)\n",
                n_multi_probe_genes))
  }

  return(expr_annotated)
}

# =====================================================
# 3. limma 差异分析
# =====================================================

#' 使用 limma 进行芯片数据差异分析
#'
#' @param expr_matrix 基因表达矩阵
#' @param ctrl_samples 对照组样本名
#' @param trt_samples 处理组样本名
#' @param pvalue_threshold P值阈值
#' @param logfc_threshold log2FC阈值
#' @param pval_type P值类型："adj.P.Val" 或 "P.Value"
#' @return 差异分析结果
run_limma_analysis <- function(expr_matrix, ctrl_samples, trt_samples,
                               pvalue_threshold = 0.05,
                               logfc_threshold = 1,
                               pval_type = "adj.P.Val") {
  cat("🧬 开始 limma 差异分析...\n")

  # 检查样本
  if (length(ctrl_samples) == 0 || length(trt_samples) == 0) {
    stop("❌ 对照组和处理组样本数不能为0")
  }

  # 重新排序矩阵（按照 Control -> Treatment）
  sample_order <- c(ctrl_samples, trt_samples)

  # 检查样本是否存在
  missing_samples <- setdiff(sample_order, colnames(expr_matrix))
  if (length(missing_samples) > 0) {
    stop(sprintf("❌ 以下样本不存在于表达矩阵中: %s",
                 paste(missing_samples, collapse = ", ")))
  }

  expr_ordered <- expr_matrix[, sample_order, drop = FALSE]

  # 🆕 检查并移除非数值列（如ProbeID和Gene）
  # limma需要纯数值表达矩阵
  if ("ProbeID" %in% colnames(expr_ordered)) {
    cat("📋 移除ProbeID列（非数值）\n")
    expr_ordered <- expr_ordered[, colnames(expr_ordered) != "ProbeID", drop = FALSE]
  }

  if ("Gene" %in% colnames(expr_ordered)) {
    cat("📋 移除Gene列（非数值）\n")
    expr_ordered <- expr_ordered[, colnames(expr_ordered) != "Gene", drop = FALSE]
  }

  # 确保所有列都是数值
  expr_ordered <- as.matrix(expr_ordered)
  storage.mode(expr_ordered) <- "numeric"

  cat(sprintf("📊 最终分析矩阵: %d 基因 × %d 样本\n",
              nrow(expr_ordered), ncol(expr_ordered)))

  # 创建分组因子
  group <- factor(c(rep("Control", length(ctrl_samples)),
                    rep("Treatment", length(trt_samples))))

  cat(sprintf("📊 样本分组: Control=%d, Treatment=%d\n",
              length(ctrl_samples), length(trt_samples)))

  # 加载 limma
  library(limma)

  # 设计矩阵
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  cat("📊 设计矩阵:\n")
  print(design)

  # 线性模型拟合
  fit <- lmFit(expr_ordered, design)

  # 设置对比
  contrast.matrix <- makeContrasts(Treatment-Control, levels=design)

  cat("📊 对比矩阵:\n")
  print(contrast.matrix)

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  # 提取结果
  results <- topTable(fit2,
                      number = Inf,
                      adjust.method = "BH",
                      sort.by = "P")

  # 🔧 修复：正确设置ID和SYMBOL列
  # rownames(results) 是基因符号（如"FAM174B", "TP53"）
  # expr_matrix现在包含EntrezID列（如果有），用于ID列
  results$SYMBOL <- rownames(results)       # SYMBOL = 基因符号（行名）

  # ID列使用EntrezID（如果存在），否则使用基因符号
  if ("EntrezID" %in% colnames(expr_matrix)) {
    # expr_matrix的行名对应results的行名
    results$ID <- expr_matrix[rownames(results), "EntrezID"]  # ID = Entrez Gene ID
    cat("✅ ID列使用Entrez Gene ID，SYMBOL列使用基因符号\n")
  } else {
    results$ID <- rownames(results)  # 没有EntrezID时，ID也用基因符号
    cat("✅ ID列和SYMBOL列都使用基因符号（无EntrezID）\n")
  }

  # 计算统计信息 - 根据用户选择的P值类型
  n_total <- nrow(results)

  # 🔧 使用用户选择的P值类型
  pval_col <- if (pval_type == "adj.P.Val") "adj.P.Val" else "P.Value"
  pval_values <- results[[pval_col]]

  n_significant <- sum(pval_values < pvalue_threshold, na.rm = TRUE)
  n_up <- sum(results$logFC > logfc_threshold & pval_values < pvalue_threshold,
              na.rm = TRUE)
  n_down <- sum(results$logFC < -logfc_threshold & pval_values < pvalue_threshold,
                na.rm = TRUE)

  cat(sprintf("✅ limma 分析完成: %d 个基因\n", n_total))
  cat(sprintf("   使用P值类型: %s\n", pval_type))
  cat(sprintf("   显著差异基因 (%s < %.3f): %d (%.1f%%)\n",
              pval_type, pvalue_threshold, n_significant, n_significant/n_total*100))
  cat(sprintf("   上调基因 (log2FC > %.2f): %d\n", logfc_threshold, n_up))
  cat(sprintf("   下调基因 (log2FC < %.2f): %d\n", logfc_threshold, n_down))

  return(list(
    results = results,
    n_total = n_total,
    n_significant = n_significant,
    n_up = n_up,
    n_down = n_down,
    design = design,
    fit = fit2
  ))
}

# =====================================================
# 4. 探针表达量聚合
# =====================================================

#' 聚合探针表达量到基因水平
#'
#' @param probe_matrix 探针表达矩阵（行为探针，列为样本）
#' @param probe_mapping 探针-基因映射数据框（从 parse_soft_platform 返回）
#' @return 基因表达矩阵（行为基因，列为样本）
aggregate_probe_expression <- function(probe_matrix, probe_mapping) {
  cat("🔄 开始聚合探针表达量...\n")

  # 提取映射关系
  probe_ids <- probe_mapping$probe_id
  gene_symbols <- probe_mapping$gene_symbol

  # 创建映射向量（探针ID -> 基因符号）
  names(gene_symbols) <- probe_ids

  # 过滤掉未映射的探针
  common_probes <- intersect(rownames(probe_matrix), probe_ids)

  if (length(common_probes) == 0) {
    cat("⚠️ 未找到匹配的探针ID\n")
    return(NULL)
  }

  cat(sprintf("✅ 匹配探针: %d / %d\n", length(common_probes), nrow(probe_matrix)))

  # 子集表达矩阵和映射
  expr_subset <- probe_matrix[common_probes, , drop = FALSE]
  gene_symbols_subset <- gene_symbols[common_probes]

  # 移除NA和空字符串
  valid_mask <- !is.na(gene_symbols_subset) & gene_symbols_subset != ""
  expr_subset <- expr_subset[valid_mask, , drop = FALSE]
  gene_symbols_subset <- gene_symbols_subset[valid_mask]

  # 聚合方法：对于每个基因，选择表达量最高的探针
  cat("📊 聚合策略: 选择最高表达探针\n")

  # 获取唯一基因
  unique_genes <- unique(gene_symbols_subset)
  cat(sprintf("📊 唯一基因数: %d\n", length(unique_genes)))

  # 为每个基因选择表达量最高的探针
  gene_expr_list <- lapply(unique_genes, function(gene) {
    # 找到该基因的所有探针
    gene_probes <- which(gene_symbols_subset == gene)

    if (length(gene_probes) == 1) {
      # 只有一个探针，直接使用
      return(expr_subset[gene_probes, , drop = FALSE])
    } else {
      # 多个探针，选择平均表达量最高的
      avg_expr <- rowMeans(expr_subset[gene_probes, , drop = FALSE])
      best_probe <- gene_probes[which.max(avg_expr)]
      return(expr_subset[best_probe, , drop = FALSE])
    }
  })

  # 合并为基因表达矩阵
  gene_expr_matrix <- do.call(rbind, gene_expr_list)
  rownames(gene_expr_matrix) <- unique_genes

  cat(sprintf("✅ 聚合完成: %d 个基因\n", nrow(gene_expr_matrix)))

  return(gene_expr_matrix)
}

# =====================================================
# 5. 智能样本分组系统（通用版本）
# =====================================================

#' 从 Series Matrix 元数据中自动检测分组
#'
#' @param sample_names 样本名向量
#' @param sample_descriptions 样本描述向量（可选）
#' @param sample_titles 样本标题向量（可选）
#' @return list 包含分组建议和分组方法
detect_chip_groups_auto <- function(sample_names,
                                     sample_descriptions = NULL,
                                     sample_titles = NULL) {
  cat("🔍 开始自动检测分组模式...\n")
  cat(sprintf("📊 总样本数: %d\n", length(sample_names)))

  # 定义常见的分组模式
  group_patterns <- list(
    # 模式1: 时间序列 (before/after, baseline/followup)
    list(
      name = "时间序列",
      ctrl_keywords = c("before", "baseline", "time0", "initial", "visit1"),
      trt_keywords = c("after", "post", "follow", "final", "visit2", "visit3")
    ),

    # 模式2: 处理对照 (control/treated)
    list(
      name = "处理对照",
      ctrl_keywords = c("control", "ctrl", "untreated", "vehicle", "placebo"),
      trt_keywords = c("treatment", "treated", "drug", "compound", "stimulated")
    ),

    # 模式3: 疾病对照 (normal/disease)
    list(
      name = "疾病对照",
      ctrl_keywords = c("normal", "healthy", "control", "wild"),
      trt_keywords = c("disease", "patient", "cancer", "tumor", "sick")
    ),

    # 模式4: 基因型 (wildtype/mutant)
    list(
      name = "基因型",
      ctrl_keywords = c("wild", "wt", "wildtype", "normal", "control"),
      trt_keywords = c("mutant", "mut", "knockout", "ko", "transgenic", "tg")
    ),

    # 模式5: 剂量反应 (control/dose)
    list(
      name = "剂量反应",
      ctrl_keywords = c("dose0", "dose_0", "control", "vehicle", "untreated"),
      trt_keywords = c("dose", "treatment", "low", "medium", "high")
    ),

    # 模式6: 激活/抑制
    list(
      name = "激活抑制",
      ctrl_keywords = c("inactive", "unstimulated", "resting", "control"),
      trt_keywords = c("active", "stimulated", "induced", "activated")
    )
  )

  # 方法1: 从 sample_descriptions 检测
  if (!is.null(sample_descriptions) && length(sample_descriptions) > 0) {
    cat("📋 尝试从 Sample_description 检测分组...\n")

    for (pattern in group_patterns) {
      ctrl_match <- sapply(sample_descriptions, function(d) {
        any(sapply(pattern$ctrl_keywords, function(kw) {
          grepl(kw, d, ignore.case = TRUE)
        }))
      })

      trt_match <- sapply(sample_descriptions, function(d) {
        any(sapply(pattern$trt_keywords, function(kw) {
          grepl(kw, d, ignore.case = TRUE)
        }))
      })

      # 检查是否找到匹配
      if (sum(ctrl_match) > 0 && sum(trt_match) > 0) {
        ctrl_idx <- which(ctrl_match)
        trt_idx <- which(trt_match)

        cat(sprintf("✅ 检测到 '%s' 分组模式 (from description)\n", pattern$name))
        cat(sprintf("   对照组: %d 个样本 (%s)\n",
                    length(ctrl_idx),
                    paste(sample_names[ctrl_idx], collapse = ", ")))
        cat(sprintf("   处理组: %d 个样本 (%s)\n",
                    length(trt_idx),
                    paste(sample_names[trt_idx], collapse = ", ")))

        return(list(
          pattern_name = pattern$name,
          method = "auto_description",
          ctrl_samples = sample_names[ctrl_idx],
          trt_samples = sample_names[trt_idx],
          ctrl_indices = ctrl_idx,
          trt_indices = trt_idx,
          confidence = "high",
          source = "description"
        ))
      }
    }
  }

  # 方法2: 从 sample_titles 检测
  if (!is.null(sample_titles) && length(sample_titles) > 0) {
    cat("📋 尝试从 Sample_title 检测分组...\n")

    for (pattern in group_patterns) {
      ctrl_match <- sapply(sample_titles, function(t) {
        any(sapply(pattern$ctrl_keywords, function(kw) {
          grepl(kw, t, ignore.case = TRUE)
        }))
      })

      trt_match <- sapply(sample_titles, function(t) {
        any(sapply(pattern$trt_keywords, function(kw) {
          grepl(kw, t, ignore.case = TRUE)
        }))
      })

      if (sum(ctrl_match) > 0 && sum(trt_match) > 0) {
        ctrl_idx <- which(ctrl_match)
        trt_idx <- which(trt_match)

        cat(sprintf("✅ 检测到 '%s' 分组模式 (from title)\n", pattern$name))
        cat(sprintf("   对照组: %d 个样本\n", length(ctrl_idx)))
        cat(sprintf("   处理组: %d 个样本\n", length(trt_idx)))

        return(list(
          pattern_name = pattern$name,
          method = "auto_title",
          ctrl_samples = sample_names[ctrl_idx],
          trt_samples = sample_names[trt_idx],
          ctrl_indices = ctrl_idx,
          trt_indices = trt_idx,
          confidence = "medium",
          source = "title"
        ))
      }
    }
  }

  # 方法3: 从样本名本身检测
  cat("📋 尝试从样本名检测分组...\n")

  # 简化样本名（移除GSM前缀）
  simplified_names <- gsub("^GSM\\d+_", "", sample_names)

  for (pattern in group_patterns) {
    ctrl_match <- sapply(simplified_names, function(n) {
      any(sapply(pattern$ctrl_keywords, function(kw) {
        grepl(kw, n, ignore.case = TRUE)
      }))
    })

    trt_match <- sapply(simplified_names, function(n) {
      any(sapply(pattern$trt_keywords, function(kw) {
        grepl(kw, n, ignore.case = TRUE)
      }))
    })

    if (sum(ctrl_match) > 0 && sum(trt_match) > 0) {
      ctrl_idx <- which(ctrl_match)
      trt_idx <- which(trt_match)

      cat(sprintf("✅ 检测到 '%s' 分组模式 (from name)\n", pattern$name))
      cat(sprintf("   对照组: %d 个样本\n", length(ctrl_idx)))
      cat(sprintf("   处理组: %d 个样本\n", length(trt_idx)))

      return(list(
        pattern_name = pattern$name,
        method = "auto_name",
        ctrl_samples = sample_names[ctrl_idx],
        trt_samples = sample_names[trt_idx],
        ctrl_indices = ctrl_idx,
        trt_indices = trt_idx,
        confidence = "medium",
        source = "name"
      ))
    }
  }

  # 所有方法都失败
  cat("⚠️ 未能自动检测到分组模式，请手动设置\n")

  return(list(
    pattern_name = NULL,
    method = "manual",
    ctrl_samples = NULL,
    trt_samples = NULL,
    ctrl_indices = NULL,
    trt_indices = NULL,
    confidence = NULL,
    source = NULL
  ))
}

#' 检测配对关系
#'
#' @param sample_names 样本名
#' @param metadata 元数据
#' @return list 配对关系或NULL
detect_pairing <- function(sample_names, metadata) {
  # 尝试从样本名中提取配对信息
  # 例如: Patient1_before, Patient1_after

  # 提取共同的配对标识符
  # 方法：移除已知的对照组/处理组关键词后，看剩余部分是否匹配

  # 简化样本名
  simplified <- gsub("^GSM\\d+_", "", sample_names)

  # 定义配对关键词
  pairing_keywords <- list(
    before = "after",
    baseline = "follow",
    control = "treated",
    time0 = "time1",
    visit1 = "visit2",
    pre = "post"
  )

  # 尝试检测配对模式
  for (kw1 in names(pairing_keywords)) {
    kw2 <- pairing_keywords[[kw1]]

    # 检查是否有样本包含这些关键词
    has_kw1 <- grepl(kw1, simplified, ignore.case = TRUE)
    has_kw2 <- grepl(kw2, simplified, ignore.case = TRUE)

    if (sum(has_kw1) > 0 && sum(has_kw2) > 0) {
      # 提取配对标识符
      ids_kw1 <- gsub(kw1, "", simplified[has_kw1], ignore.case = TRUE)
      ids_kw2 <- gsub(kw2, "", simplified[has_kw2], ignore.case = TRUE)

      # 找到共同的ID
      common_ids <- intersect(ids_kw1, ids_kw2)

      if (length(common_ids) > 0) {
        cat(sprintf("💡 检测到配对设计: %d 对样本\n", length(common_ids)))

        return(list(
          is_paired = TRUE,
          n_pairs = length(common_ids),
          pairing_pattern = sprintf("%s/%s", kw1, kw2)
        ))
      }
    }
  }

  return(list(is_paired = FALSE))
}

# =====================================================
# 5. 转换为标准格式（与现有模块兼容）
# =====================================================

#' 将芯片差异分析结果转换为标准格式
#'
#' @param limma_results limma 分析结果
#' @param expr_matrix 表达矩阵
#' @param ctrl_samples 对照组样本
#' @param trt_samples 处理组样本
#' @return 标准格式的差异分析结果
format_chip_results_for_pipeline <- function(limma_results, expr_matrix,
                                             ctrl_samples, trt_samples) {
  # 🔧 修复：正确使用ID（EntrezID）和SYMBOL（基因符号）列
  limma_res <- limma_results$results

  deg_df <- data.frame(
    ID = limma_res$ID,                          # Entrez Gene ID（如果有的话）
    SYMBOL = limma_res$SYMBOL,                  # 基因符号
    log2FoldChange = limma_res$logFC,           # log2倍数变化
    pvalue = limma_res$P.Value,                 # 原始p值
    padj = limma_res$adj.P.Val,                 # BH校正p值
    baseMean = limma_res$AveExpr,               # 平均表达
    t = limma_res$t,                            # t统计量（用于TF活性分析）
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # ENTREZID列与ID列相同（都是Entrez Gene ID）
  deg_df$ENTREZID <- deg_df$ID

  return(list(
    deg_df = deg_df,
    background_genes = rownames(expr_matrix),   # 背景基因（用于富集分析）
    expr_matrix = expr_matrix,                  # 完整表达矩阵（用于 AUCell/GSVA）
    ctrl_samples = ctrl_samples,
    trt_samples = trt_samples,
    method = "limma",
    n_significant = limma_results$n_significant,
    n_up = limma_results$n_up,
    n_down = limma_results$n_down
  ))
}

# =====================================================
# 6. Shiny Server 函数
# =====================================================

#' 芯片数据分析模块 UI 生成函数
#'
#' @return UI 元素
chip_analysis_ui <- function() {
  tagList(
    fluidRow(
      column(12,
        div(
          class = "info-box",
          style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                 color: white; padding: 25px; border-radius: 15px; margin-bottom: 25px;",
          h4("🧬 芯片数据分析模块", style = "margin-top: 0; color: white;"),
          p("支持 GEO Series Matrix 格式的芯片数据差异分析，自动探针注释，无缝集成下游富集分析。",
            style = "color: rgba(255,255,255,0.9); margin-bottom: 0;")
        )
      )
    ),

    # 🎨 使用可折叠面板组织所有步骤
    tags$div(
      id = "chip_analysis_accordion",

      # ===== 面板1: 数据上传 =====
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          style = "cursor: pointer; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px;",
          `data-toggle` = "collapse",
          `data-target` = "#panel_upload",
          tags$h4(
            class = "panel-title",
            style = "margin: 0;",
            tags$span(icon("upload"), " 📁 步骤1: 数据上传与预览")
          )
        ),
        tags$div(
          id = "panel_upload",
          class = "panel-collapse collapse in",  # in = 默认展开

        wellPanel(
          # 文件上传
          fluidRow(
            column(6,
              h5("📄 上传数据文件", style = "color: #007AFF;"),
              fileInput("chip_series_matrix",
                       "GEO Series Matrix 文件",
                       accept = c(".txt", ".matrix.txt", "text/plain"),
                       placeholder = "选择文件..."),
              helpText("GEO 数据库下载的 Series Matrix 文件（通常包含样本表达矩阵）")
            ),
            column(6,
              fileInput("chip_soft_platform",
                       "SOFT 平台注释文件 (可选)",
                       accept = c(".txt", ".soft", "annot.txt", "text/plain"),
                       placeholder = "选择文件..."),
              helpText("用于探针注释的 GPL 平台文件。如不上传，系统将尝试自动注释。")
            )
          ),

          tags$hr(style = "border-color: #dee2e6;"),

          # 数据预览
          h5("📊 数据文件预览", style="color: #007AFF;"),

          fluidRow(
            column(6,
              h6("Series Matrix 文件（前5行）", style="color: #666;"),
              DTOutput("chip_series_matrix_preview")
            ),
            column(6,
              h6("SOFT 文件（前10行）", style="color: #666;"),
              DTOutput("chip_soft_raw_preview")
            )
          )
        ),
        tags$div(
          class = "panel-body",
          style = "padding: 15px;"
        )
      )
    ),

    # ===== 面板2: 探针注释配置 =====
    tags$div(
      class = "panel panel-default",
      tags$div(
        class = "panel-heading",
        style = "cursor: pointer; background: linear-gradient(135deg, #9C27B0 0%, #7B1FA2 100%); color: white; padding: 15px;",
        `data-toggle` = "collapse",
        `data-target` = "#panel_annotation",
        tags$h4(
          class = "panel-title",
          style = "margin: 0;",
          tags$span(icon("cogs"), " 🧬 步骤2: 探针注释与数据合并")
        )
      ),
      tags$div(
        id = "panel_annotation",
        class = "panel-collapse collapse",  # 默认折叠

        wellPanel(
          # SOFT文件列名清单
          uiOutput("chip_soft_columns_list_ui"),

          tags$hr(style = "border-color: #dee2e6;"),

          # 探针注释配置
          h5("📋 探针注释配置", style = "color: #9C27B0;"),
          uiOutput("chip_soft_column_selection_panel"),

          tags$hr(style = "border-color: #dee2e6;"),

          # 合并操作按钮
          fluidRow(
            column(6,
              actionButton("chip_preview_merge", "👁️ 预览合并结果",
                          class = "btn-info", style = "width: 100%;")
            ),
            column(6,
              actionButton("chip_apply_merge", "✅ 应用配置并生成最终矩阵",
                          class = "btn-success", style = "width: 100%;")
            )
          ),

          # 合并预览
          conditionalPanel(
            condition = "input.chip_preview_merge",
            wellPanel(
              style = "background: #e8f5e9; border: 2px solid #4caf50;",
              h5("👁️ 合并结果预览（前5行）", style = "color: #2e7d32;"),
              DTOutput("chip_merge_preview_table")
            )
          ),

          # 最终矩阵显示
          conditionalPanel(
            condition = "input.chip_apply_merge",
            wellPanel(
              style = "background: linear-gradient(135deg, #e8f5e9 0%, #c8e6c9 100%); border: 2px solid #4caf50;",
              h5("✅ 最终表达矩阵", style = "color: #2e7d32; font-size: 18px; font-weight: bold;"),
              uiOutput("chip_final_matrix_ui"),
              br(),
              helpText("💡 此矩阵可直接用于后续的差异分析。已将探针ID和基因符号合并到表达数据中。", style = "color: #2e7d33;")
            )
          )
        )
      ),

      # ===== 面板3: 数据预处理与探针去重 =====
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          style = "cursor: pointer; background: linear-gradient(135deg, #ff9800 0%, #f57c00 100%); color: white; padding: 15px;",
          `data-toggle` = "collapse",
          `data-target` = "#panel_preprocess",
          tags$h4(
            class = "panel-title",
            style = "margin: 0;",
            tags$span(icon("sliders-h"), " 🔧 步骤3: 数据预处理与探针去重")
          )
        ),
        tags$div(
          id = "panel_preprocess",
          class = "panel-collapse collapse in",  # 默认展开

        wellPanel(
          # 预处理
          h5("📊 数据预处理（log2转换 + 标准化）", style = "color: #ff9800; font-weight: bold;"),
          fluidRow(
            column(6,
              checkboxInput("chip_auto_log2", "自动判断并执行log2转换", value = TRUE),
              checkboxInput("chip_normalize_data", "执行limma标准化（normalizeBetweenArrays）", value = TRUE)
            ),
            column(6,
              actionButton("chip_preprocess_data", "⚙️ 执行预处理",
                          class = "btn-warning", style = "width: 100%;")
            )
          ),

          # 预处理结果
          conditionalPanel(
            condition = "input.chip_preprocess_data",
            uiOutput("chip_preprocess_result_ui")
          ),

          tags$hr(style = "border-color: #ffc107;"),

          # 批次矫正
          h5("🎛️ 批次效应矫正（可选）", style = "color: #E91E63; font-weight: bold;"),
          fluidRow(
            column(6,
              selectInput("chip_batch_method", "批次矫正方法",
                          choices = c("无" = "none", "ComBat (sva)" = "combat",
                                    "ComBat (limma)" = "limma", "SVA" = "sva"),
                          selected = "none")
            ),
            column(6,
              actionButton("chip_apply_batch_correct", "🎛️ 执行批次矫正",
                          class = "btn-danger", style = "width: 100%;")
            )
          ),

          # 批次矫正结果
          conditionalPanel(
            condition = "input.chip_apply_batch_correct",
            uiOutput("chip_batch_correct_result_ui")
          ),

          tags$hr(style = "border-color: #ffc107;"),

          # 探针去重
          h5("✂️ 探针去重（保留表达量最高的探针）", style = "color: #ff9800; font-weight: bold;"),
          helpText("当一个基因对应多个探针时，保留表达量最高的探针。这将生成基因级别的表达矩阵。"),

          fluidRow(
            column(6,
              h6("去重前统计：", style = "color: #666;"),
              uiOutput("chip_before_dedupe_stats")
            ),
            column(6,
              h6("去重后统计：", style = "color: #666;"),
              uiOutput("chip_after_dedupe_stats")
            )
          ),

          actionButton("chip_dedupe_probes", "✂️ 执行探针去重",
                      class = "btn-warning btn-lg", style = "width: 100%;"),

          # 去重结果
          conditionalPanel(
            condition = "input.chip_dedupe_probes",
            wellPanel(
              h6("✅ 去重完成", style = "color: #28a745;"),
              uiOutput("chip_dedupe_result_ui")
            )
          ),

          tags$hr(style = "border-color: #ffc107;"),

          # 生成标准格式数据
          h5("💾 生成标准格式数据", style = "color: #ff9800; font-weight: bold;"),
          helpText("将处理后的表达矩阵转换成标准格式，可直接用于后续的差异分析、KEGG、GO等模块。"),

          fluidRow(
            column(12,
              div(
                style = "background: #d4edda; padding: 15px; border-radius: 8px; border: 1px solid #c3e6cb;",
                h6("💡 即可对接现有分析模块：", style = "color: #155724;"),
                tags$ul(style="margin: 10px 0; padding-left: 20px;",
                  tags$li("差异分析（使用现有的差异分析模块）"),
                  tags$li("KEGG富集分析（使用现有的KEGG模块）"),
                  tags$li("GO富集分析（使用现有的GO模块）"),
                  tags$li("GSEA分析（使用现有的GSEA模块）")
                )
              )
            )
          ),

          actionButton("chip_generate_standard_data", "🚀 生成标准格式数据",
                      class = "btn-success btn-lg", style = "width: 100%; font-size: 16px;")
          )
        )
    ),

      # ===== 面板4: 差异分析 =====
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          style = "cursor: pointer; background: linear-gradient(135deg, #34C759 0%, #2e7d32 100%); color: white; padding: 15px;",
          `data-toggle` = "collapse",
          `data-target` = "#panel_diff_analysis",
          tags$h4(
            class = "panel-title",
            style = "margin: 0;",
            tags$span(icon("chart-bar"), " 🧬 步骤4: 差异分析")
          )
        ),
        tags$div(
          id = "panel_diff_analysis",
          class = "panel-collapse collapse",  # 默认折叠
          style = "padding: 15px;",

          wellPanel(
            # 样本分组
            uiOutput("chip_grouping_ui"),

            tags$hr(style = "border-color: #34C759;"),

            # 差异分析参数
            h5("🔬 差异分析参数", style = "color: #34C759; font-weight: bold;"),
            fluidRow(
              column(4,
                sliderInput("chip_logfc_threshold",
                            "log2FoldChange 阈值:",
                            min = 0, max = 5, value = 1, step = 0.1)
              ),
              column(4,
                selectInput("chip_pval_type",
                            "显著性指标:",
                            choices = c("校正P值 (adj.P.Val)" = "adj.P.Val",
                                        "原始P值 (P.Value)" = "P.Value"),
                            selected = "adj.P.Val")
              ),
              column(4,
                sliderInput("chip_pvalue_threshold",
                            "P值 阈值:",
                            min = 0.001, max = 0.1, value = 0.05, step = 0.001)
              )
            ),

            fluidRow(
              column(12,
                checkboxInput("chip_paired_analysis",
                             "配对样本分析（如果适用）",
                             value = FALSE)
              )
            ),

            tags$hr(style = "border-color: #34C759;"),

            actionButton("run_chip_analysis", "🚀 运行差异分析",
                        class = "btn-primary btn-lg",
                        style = "width: 100%; margin-top: 15px;")
          )
        )
      )
    ),

      # ===== 面板5: 分析结果 =====
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          style = "cursor: pointer; background: linear-gradient(135deg, #007AFF 0%, #0051D5 100%); color: white; padding: 15px;",
          `data-toggle` = "collapse",
          `data-target` = "#panel_results",
          tags$h4(
            class = "panel-title",
            style = "margin: 0;",
            tags$span(icon("table"), " 📊 步骤5: 分析结果")
          )
        ),
        tags$div(
          id = "panel_results",
          class = "panel-collapse collapse",  # 默认折叠
          style = "padding: 15px;",

          wellPanel(
            # 结果统计和表格（chip_results_ui 已包含表格）
            uiOutput("chip_results_ui"),

            tags$hr(style = "border-color: #007AFF;"),

            # 下载按钮
            fluidRow(
              column(12,
                downloadButton("download_chip_results", "📥 下载结果", class = "btn-success")
              )
            )
          )
        )
      )
    )
  )
}

#' 芯片数据分析模块 Server 函数
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param deg_results 差异分析结果容器（用于更新）
chip_analysis_server <- function(input, output, session, deg_results) {
  cat("✅ 芯片分析模块已启动\n")

  # 渲染 UI
  output$chip_analysis_ui_output <- renderUI({
    chip_analysis_ui()
  })

  # 存储芯片数据
  chip_data <- reactiveValues(
    series_matrix = NULL,
    soft_platform = NULL,
    probe_mapping = NULL,
    expr_matrix = NULL,
    metadata = NULL,
    group_info = NULL,
    manual_ctrl_samples = NULL,
    manual_trt_samples = NULL
  )

  # 解析对照组样本
  observeEvent(input$chip_parse_ctrl, {
    req(input$chip_paste_ctrl)
    req(chip_data$series_matrix)

    pasted_text <- input$chip_paste_ctrl

    # 显示进度
    showNotification("正在解析对照组样本...", type = "message")

    # 解析样本列表
    samples <- parse_sample_list(pasted_text, chip_data$series_matrix)

    if (is.null(samples) || length(samples) == 0) {
      showNotification("解析失败：未找到有效样本", type = "error")
      return(NULL)
    }

    # 保存对照组
    chip_data$manual_ctrl_samples <- samples

    showNotification(
      sprintf("✅ 成功解析对照组: %d 个样本", length(samples)),
      type = "message"
    )
  })

  # 解析处理组样本
  observeEvent(input$chip_parse_trt, {
    req(input$chip_paste_trt)
    req(chip_data$series_matrix)

    pasted_text <- input$chip_paste_trt

    # 显示进度
    showNotification("正在解析处理组样本...", type = "message")

    # 解析样本列表
    samples <- parse_sample_list(pasted_text, chip_data$series_matrix)

    if (is.null(samples) || length(samples) == 0) {
      showNotification("解析失败：未找到有效样本", type = "error")
      return(NULL)
    }

    # 保存处理组
    chip_data$manual_trt_samples <- samples

    showNotification(
      sprintf("✅ 成功解析处理组: %d 个样本", length(samples)),
      type = "message"
    )
  })

  # 清除分组设置
  observeEvent(input$chip_clear_groups, {
    chip_data$manual_ctrl_samples <- NULL
    chip_data$manual_trt_samples <- NULL
    showNotification("已清除分组设置", type = "message")
  })

  # 解析 Series Matrix 文件
  observeEvent(input$chip_series_matrix, {
    req(input$chip_series_matrix)

    file_path <- input$chip_series_matrix$datapath

    # 显示进度
    showNotification("正在解析 GEO Series Matrix 文件...", type = "message")

    # 解析文件
    result <- parse_geo_series_matrix(file_path)

    if (!result$success) {
      showNotification(result$error, type = "error")
      return(NULL)
    }

    # 保存数据
    chip_data$series_matrix <- result$matrix
    chip_data$metadata <- result$metadata
    chip_data$expr_matrix <- result$matrix

    showNotification(
      sprintf("✅ 成功解析: %d 探针 × %d 样本",
              result$n_probes, result$n_samples),
      type = "message"
    )

    # 尝试自动检测分组
    if (!is.null(result$metadata)) {
      group_info <- detect_chip_groups_auto(
        sample_names = result$sample_names,
        sample_descriptions = result$metadata$sample_descriptions,
        sample_titles = result$metadata$sample_titles
      )
      chip_data$group_info <- group_info
    }
  }, ignoreNULL = TRUE)

  # Series Matrix 文件预览
  output$chip_series_matrix_preview <- renderDT({
    req(chip_data$series_matrix)

    # 提取前5行
    preview_matrix <- head(chip_data$series_matrix, 5)

    # 转换为数据框（保留行名）
    preview_df <- as.data.frame(preview_matrix)

    datatable(
      preview_df,
      options = list(
        dom = 't',
        paging = FALSE,
        scrollX = TRUE,
        columnDefs = list(list(
          className = 'dt-center',
          targets = "_all"
        ))
      ),
      rownames = TRUE,  # 显示行名（探针ID）
      filter = 'none'
    ) %>%
      formatStyle(columns = 1:min(5, ncol(preview_df)), fontSize = '85%')
  })

  # SOFT 文件预览UI
  output$chip_soft_preview_ui <- renderUI({
    # 只要有文件上传记录就显示预览
    if (!is.null(input$chip_soft_platform)) {
      cat(sprintf("🔍 渲染SOFT预览UI (文件已上传)\n"))

      if (!is.null(chip_data$soft_platform)) {
        cat(sprintf("  数据可用: %d rows x %d cols\n",
                    nrow(chip_data$soft_platform), ncol(chip_data$soft_platform)))
      } else {
        cat("  ⚠️ 数据尚未加载到chip_data\n")
      }

      tagList(
        h5("📄 SOFT 平台注释文件（前10行）", style = "color: #FF9800;"),
        helpText("这是SOFT注释文件的真实数据。可以看到ID列和多个可能的基因列。"),
        DTOutput("chip_soft_raw_preview")
      )
    } else {
      cat("⚠️ SOFT文件未上传\n")
      div(
        class = "alert alert-info",
        h5("📄 请先上传SOFT平台注释文件"),
        p("上传后将在此处显示文件预览。")
      )
    }
  })

  # 🆕 SOFT文件列名列表展示
  output$chip_soft_columns_list_ui <- renderUI({
    req(chip_data$soft_platform)

    soft_cols <- colnames(chip_data$soft_platform)

    tagList(
      h6("📌 所有列名（共", span(style = "color: #FF9500;", length(soft_cols)), "列）:", style = "color: #333;"),

      # 以表格形式展示所有列名
      div(
        style = "background: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #FF9800;",

        # 每行显示3列
        div(style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px;",
          lapply(seq_along(soft_cols), function(i) {
            col_name <- soft_cols[i]
            div(
              style = "background: white; padding: 8px; border-radius: 4px; border: 1px solid #FFB74D;",
              tags$span(style = "color: #E65100; font-weight: bold; font-family: monospace;",
                       sprintf("%d. %s", i, col_name))
            )
          })
        )
      ),

      br(),

      # 重要提示
      div(
        style = "background: #e3f2fd; padding: 12px; border-radius: 5px; border-left: 4px solid #2196F3;",
        h6("💡 如何选择基因列？", style = "color: #1976D2; margin-top: 0;"),
        tags$ul(style = "padding-left: 20px; margin: 5px 0;",
          tags$li("查看下方的列内容示例，了解每列包含什么数据"),
          tags$li("基因符号列通常包含：TP53, EGFR, BRCA1, MYC 等基因名称"),
          tags$li("ID列通常包含：数字或探针标识符（如 1553601_at）"),
          tags$li("点击下方的'查看列内容'按钮查看每列的前5行数据")
        )
      )
    )
  })

  # 🆕 显示各列内容示例
  output$chip_soft_column_examples_ui <- renderUI({
    req(chip_data$soft_platform)

    soft_cols <- colnames(chip_data$soft_platform)

    # 为每一列生成示例展示
    tagList(
      h6("📊 各列内容示例（前3行）", style = "color: #333;"),

      div(style = "max-height: 400px; overflow-y: auto;"),
      lapply(soft_cols, function(col_name) {
        # 获取该列的前3个非空值
        col_data <- chip_data$soft_platform[[col_name]]
        col_data <- col_data[!is.na(col_data) & col_data != ""]
        examples <- head(col_data, 3)

        wellPanel(
          style = "padding: 10px; margin-bottom: 10px;",
          h7(style = "color: #FF9800; font-weight: bold; margin-bottom: 5px;",
             sprintf("🔹 %s", col_name)),
          div(
            style = "background: #f5f5f5; padding: 8px; border-radius: 4px; font-family: monospace; font-size: 11px;",
            for (i in seq_along(examples)) {
              tags$div(
                sprintf("  %d. %s", i, as.character(examples[i])),
                style = i < length(examples) ? "margin-bottom: 5px;" : ""
              )
            }
          )
        )
      })
    )
  })

  # 解析 SOFT 平台文件
  observeEvent(input$chip_soft_platform, {
    req(input$chip_soft_platform)

    file_path <- input$chip_soft_platform$datapath

    showNotification("正在解析 SOFT 平台注释文件...", type = "message")

    result <- parse_platform_annotation(file_path, "\t")  # 使用Tab作为分隔符

    if (!result$success) {
      showNotification(result$error, type = "error")
      return(NULL)
    }

    chip_data$soft_platform <- result$raw_table
    chip_data$probe_mapping <- NULL  # 不使用自动映射
    chip_data$gene_symbol_col <- NULL  # 不使用自动检测的列名

    cat(sprintf("💾 SOFT数据已保存: %d rows x %d cols\n",
                nrow(chip_data$soft_platform), ncol(chip_data$soft_platform)))
    cat("⚠️ 请用户手动选择ID列和基因列\n")

    showNotification(
      "✅ SOFT文件已加载，请在下方手动选择ID列和基因列",
      type = "message"
    )
  }, ignoreNULL = TRUE)

  # 重新解析SOFT文件（修改分隔符后）
  observeEvent(input$chip_reparse_soft, {
    req(input$chip_soft_platform)

    file_path <- input$chip_soft_platform$datapath

    # 获取分隔符，默认为Tab
    separator <- if (is.null(input$chip_soft_separator) || input$chip_soft_separator == "") {
      "\t"
    } else {
      input$chip_soft_separator
    }

    showNotification("正在重新解析 SOFT 平台注释文件...", type = "message")

    result <- parse_platform_annotation(file_path, separator)

    if (!result$success) {
      showNotification(result$error, type = "error")
      return(NULL)
    }

    chip_data$soft_platform <- result$raw_table
    chip_data$probe_mapping <- NULL  # 不使用自动映射
    chip_data$gene_symbol_col <- NULL  # 不使用自动检测的列名

    showNotification(
      "✅ 重新解析完成，请手动选择ID列和基因列",
      type = "message"
    )
  })

  # 数据概览
  output$chip_data_summary <- renderDT({
    req(chip_data$series_matrix)

    matrix <- chip_data$series_matrix

    # 创建摘要表
    summary_df <- data.frame(
      项目 = c("探针数", "样本数", "样本名称"),
      值 = c(
        nrow(matrix),
        ncol(matrix),
        paste(colnames(matrix), collapse = ", ")
      )
    )

    # 添加分组信息
    if (!is.null(chip_data$group_info) && !is.null(chip_data$group_info$pattern_name)) {
      summary_df <- rbind(summary_df, data.frame(
        项目 = c("检测到分组模式", "对照组样本", "处理组样本"),
        值 = c(
          chip_data$group_info$pattern_name,
          paste(chip_data$group_info$ctrl_samples, collapse = ", "),
          paste(chip_data$group_info$trt_samples, collapse = ", ")
        )
      ))
    }

    datatable(summary_df,
              options = list(dom = 't', paging = FALSE),
              rownames = FALSE)
  })

  # 注释状态显示
  output$chip_annotation_status <- renderUI({
    req(chip_data$series_matrix)

    # 检查是否加载了SOFT文件
    soft_loaded <- !is.null(chip_data$soft_platform)

    if (soft_loaded) {
      # 已加载SOFT文件但尚未配置映射
      div(
        class = "alert alert-info",
        h5("✅ SOFT文件已加载，等待配置", style = "color: #17a2b8;"),
        p(sprintf("总探针数: %d", nrow(chip_data$series_matrix))),
        p(sprintf("SOFT平台数据: %d 行 x %d 列",
                nrow(chip_data$soft_platform),
                ncol(chip_data$soft_platform))),
        p("💡 请在下方选择ID列和基因列以建立探针映射。", style = "color: #007bff; font-weight: bold;"),
        if (!is.null(chip_data$probe_mapping)) {
          p(sprintf("成功映射: %d (%.1f%%)",
                  nrow(chip_data$probe_mapping),
                  nrow(chip_data$probe_mapping) / nrow(chip_data$series_matrix) * 100))
        }
      )
    } else {
      # 未加载SOFT文件
      div(
        class = "alert alert-warning",
        h5("⚠️ 未加载探针注释文件", style = "color: #ffc107;"),
        p("将直接使用探针ID作为基因符号进行分析。"),
        p("强烈建议上传 SOFT 平台注释文件以获得准确的结果。")
      )
    }
  })

  # 🆕 SOFT文件列选择面板（使用renderUI而不是conditionalPanel）
  output$chip_soft_column_selection_panel <- renderUI({
    # 检查SOFT文件是否加载
    if (is.null(chip_data$soft_platform)) {
      return(NULL)
    }

    # 只打印一次日志，避免刷屏
    if (is.null(chip_data$panel_initialized)) {
      cat(sprintf("✅ 初始化SOFT列选择面板: %d 行 x %d 列\n",
                  nrow(chip_data$soft_platform),
                  ncol(chip_data$soft_platform)))
      chip_data$panel_initialized <- TRUE
    }

    # 使用isolate防止input变化触发重新渲染
    soft_cols <- isolate(colnames(chip_data$soft_platform))

    # 获取当前选择值（使用isolate避免建立依赖）
    current_id <- isolate({
      if (!is.null(input$chip_soft_id_col) && input$chip_soft_id_col != "") {
        input$chip_soft_id_col
      } else {
        ""
      }
    })

    current_gene <- isolate({
      if (!is.null(input$chip_soft_gene_col) && input$chip_soft_gene_col != "") {
        input$chip_soft_gene_col
      } else {
        ""
      }
    })

    wellPanel(
      style = "background: linear-gradient(135deg, #fff7e6 0%, #ffe6b3 100%); border: 2px solid #FF9800;",

      h4("📋 SOFT文件列选择", style = "color: #FF9800; margin-top: 0;"),

      helpText("💡 您可以预先浏览和选择SOFT文件的列，即使还未上传Series Matrix文件。"),

      br(),

      fluidRow(
        column(4,
          h5("选择ID列", style = "color: #9C27B0; font-weight: bold;"),
          selectInput("chip_soft_id_col",
                     "选择ID列（必须与Series Matrix的探针ID匹配）",
                     choices = c("", soft_cols),
                     selected = current_id),
          helpText("选择SOFT文件中包含探针ID的列（通常为'ID'列）")
        ),
        column(4,
          h5("选择基因列", style = "color: #9C27B0; font-weight: bold;"),
          selectInput("chip_soft_gene_col",
                     "选择基因列（包含基因符号的列）",
                     choices = c("", soft_cols),
                     selected = current_gene),
          helpText("选择包含基因符号的列（如GENE_SYMBOL, SYMBOL）")
        ),
        column(4,
          h5("选择EntrezID列（可选）", style = "color: #FF5722; font-weight: bold;"),
          selectInput("chip_soft_entrez_col",
                     "选择Entrez Gene ID列",
                     choices = c("", "自动检测", soft_cols),
                     selected = ""),
          helpText("选择包含Entrez Gene ID的列（如GENE, ENTREZID）", style = "font-size: 11px; color: #FF5722;")
        )
      ),

      # 状态提示单独的uiOutput，避免循环
      uiOutput("chip_selection_status")
    )
  })

  # 🆕 状态提示单独渲染（避免导致父级重新渲染）
  output$chip_selection_status <- renderUI({
    # 使用isolate避免触发父级renderUI
    id_col <- isolate(input$chip_soft_id_col)
    gene_col <- isolate(input$chip_soft_gene_col)
    entrez_col <- isolate(input$chip_soft_entrez_col)

    # 🔍 自动检测可能的EntrezID列
    possible_entrez_hint <- ""
    if (!is.null(chip_data$soft_platform)) {
      soft_cols <- colnames(chip_data$soft_platform)
      # 检查是否有常见的EntrezID列名
      entrez_candidates <- c("ENTREZ_GENE_ID", "ENTREZID", "EntrezID", "GeneID",
                            "GENE_ID", "ENTREZ_GENE", "GENE")
      found_cols <- intersect(entrez_candidates, soft_cols)

      if (length(found_cols) > 0) {
        possible_entrez_hint <- p(sprintf("💡 检测到可能的EntrezID列: %s",
                                          paste(found_cols, collapse = ", ")),
                                  collapse = " ")
      }
    }

    if (is.null(id_col) || id_col == "" || is.null(gene_col) || gene_col == "") {
      div(
        style = "background: #fff3cd; padding: 10px; border-radius: 5px; margin-top: 15px;",
        h6("⚠️ 未完成选择", style = "color: #856404;"),
        p("请选择ID列和基因列，然后上传Series Matrix文件并点击应用按钮。", style = "font-size: 12px;"),
        if (possible_entrez_hint != "") {
          p(possible_entrez_hint, style = "font-size: 11px; color: #FF5722; font-weight: bold;")
        }
      )
    } else {
      status_color <- "d4edda"
      status_text <- "✅ 已选择列"
      status_text_color <- "155724"

      # 检查是否选择了EntrezID列
      if (is.null(entrez_col) || entrez_col == "") {
        if (possible_entrez_hint != "") {
          status_color <- "fff3cd"
          status_text <- "⚠️ 建议选择EntrezID列"
          status_text_color <- "856404"
        }
      }

      div(
        style = sprintf("background: %s; padding: 10px; border-radius: 5px; margin-top: 15px;", status_color),
        h6(status_text, style = sprintf("color: %s;", status_text_color)),
        p(sprintf("ID列: %s | 基因列: %s", id_col, gene_col),
          style = "font-size: 12px; font-weight: bold;"),
        if (!is.null(entrez_col) && entrez_col != "" && entrez_col != "自动检测") {
          p(sprintf("EntrezID列: %s", entrez_col),
            style = "font-size: 12px; color: #FF5722; font-weight: bold;")
        },
        p("💡 请上传Series Matrix文件，然后点击'✅ 应用配置并生成最终矩阵'按钮。",
          style = "font-size: 12px;"),
        if (possible_entrez_hint != "" && (is.null(entrez_col) || entrez_col == "")) {
          p(possible_entrez_hint,
            style = "font-size: 11px; color: #FF5722; font-weight: bold; margin-top: 5px;"
          )
        }
      )
    }
  })

  # 🆕 监听列选择（只在日志中记录，不触发UI更新）
  observe({
    # 只有两个列都选择了才处理
    req(input$chip_soft_id_col)
    req(input$chip_soft_gene_col)
    req(input$chip_soft_id_col != "")
    req(input$chip_soft_gene_col != "")

    # 保存选择到chip_data
    chip_data$selected_id_col <- input$chip_soft_id_col
    chip_data$selected_gene_col <- input$chip_soft_gene_col

    # 使用isolate避免触发UI重新渲染
    isolate({
      cat(sprintf("📋 用户已选择: ID列=%s, 基因列=%s\n",
                  input$chip_soft_id_col,
                  input$chip_soft_gene_col))
    })
  })

  # 🆕 合并工作流面板renderUI（替代conditionalPanel）
  output$chip_merge_workflow_panel <- renderUI({
    # 检查两个文件是否都已上传
    has_series <- !is.null(chip_data$series_matrix)
    has_soft <- !is.null(chip_data$soft_platform)

    if (!has_series || !has_soft) {
      return(NULL)
    }

    # 检查是否已选择列
    has_id_col <- !is.null(input$chip_soft_id_col) && input$chip_soft_id_col != ""
    has_gene_col <- !is.null(input$chip_soft_gene_col) && input$chip_soft_gene_col != ""

    wellPanel(
      style = "background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%); border: 2px solid #667eea;",

      h4("🔗 探针注释与数据合并工作流", style = "color: #667eea; margin-top: 0;"),

      # 工作流说明
      div(
        style = "background: rgba(255,255,255,0.8); padding: 15px; border-radius: 8px; margin-bottom: 20px;",
        h5("📋 重要提示", style = "color: #dc3545;"),
        helpText("您已在上方黄色面板选择了列，现在可以应用配置生成最终矩阵！"),
        if (has_id_col && has_gene_col) {
          tags$div(style = "background: #d4edda; padding: 10px; border-radius: 5px; border-left: 4px solid #28a745;",
            tags$strong("✅ 当前配置："),
            tags$ul(style = "margin: 10px 0;",
              tags$li(sprintf("ID列: %s", input$chip_soft_id_col)),
              tags$li(sprintf("基因列: %s", input$chip_soft_gene_col))
            )
          )
        } else {
          tags$div(style = "background: #fff3cd; padding: 10px; border-radius: 5px; border-left: 4px solid #ffc107;",
            tags$strong("⚠️ 请先在上方黄色面板选择列！")
          )
        }
      ),

      # 应用按钮
      if (has_id_col && has_gene_col) {
        tagList(
          hr(style = "border-color: #667eea;"),
          h5("步骤5: 应用配置并生成最终矩阵", style = "color: #9C27B0; font-weight: bold;"),
          fluidRow(
            column(12,
              actionButton("chip_apply_merge", "✅ 应用配置并生成最终矩阵",
                          class = "btn-success btn-lg btn-block",
                          style = "font-size: 16px; padding: 15px;"),
              helpText("点击后将应用所有配置，生成带基因符号的表达矩阵。", style = "text-align: center;")
            )
          )
        )
      }
    )
  })

  # SOFT文件列选择UI
  output$chip_soft_columns_ui <- renderUI({
    req(chip_data$soft_platform)

    soft_cols <- colnames(chip_data$soft_platform)

    tagList(
      h5("🔍 SOFT文件列信息", style = "color: #9C27B0;"),
      p(sprintf("检测到 %d 列，请确认基因符号列：", ncol(chip_data$soft_platform))),
      fluidRow(
        column(8,
          selectInput("chip_gene_symbol_col",
                     "选择基因符号列:",
                     choices = soft_cols,
                     selected = chip_data$gene_symbol_col)
        ),
        column(4,
          br(),
          actionButton("chip_update_gene_col", "🔄 更新列选择",
                      class = "btn-primary", style = "width: 100%;")
        )
      ),
      helpText("💡 提示：系统已自动检测，如不正确可手动选择。",
               class = "text-info")
    )
  })

  # 更新基因符号列
  observeEvent(input$chip_update_gene_col, {
    req(chip_data$soft_platform)
    req(input$chip_gene_symbol_col)

    selected_col <- input$chip_gene_symbol_col

    # 重新提取映射
    probe_col <- chip_data$soft_platform[, 1]  # 假设第一列是探针ID
    gene_col <- chip_data$soft_platform[, selected_col]

    # 移除NA和空
    valid_mask <- !is.na(gene_col) & gene_col != ""
    probe_col <- probe_col[valid_mask]
    gene_col <- gene_col[valid_mask]

    # 创建映射
    mapping <- data.frame(
      probe_id = as.character(probe_col),
      gene_symbol = as.character(gene_col),
      stringsAsFactors = FALSE
    )

    chip_data$probe_mapping <- mapping
    chip_data$gene_symbol_col <- selected_col

    showNotification(
      sprintf("✅ 已更新基因列为: %s (%d 个映射)", selected_col, nrow(mapping)),
      type = "message"
    )
  })

  # SOFT文件原始数据预览
  output$chip_soft_raw_preview <- renderDT({
    req(chip_data$soft_platform)

    # 显示前10行，但限制显示的列数
    preview_df <- head(chip_data$soft_platform, 10)

    # 🔧 截断过长的文本以避免渲染问题
    truncate_text <- function(text, max_len = 100) {
      if (is.character(text) || is.factor(text)) {
        text <- as.character(text)
        text <- ifelse(nchar(text) > max_len,
                       paste0(substr(text, 1, max_len), "..."),
                       text)
      }
      return(text)
    }

    # 对所有列应用截断
    for (col in colnames(preview_df)) {
      preview_df[[col]] <- truncate_text(preview_df[[col]], max_len = 100)
    }

    # 限制显示的列数（最多显示15列，避免表格过宽）
    if (ncol(preview_df) > 15) {
      preview_df <- preview_df[, 1:15]
      cat(sprintf("⚠️ SOFT文件有%d列，仅显示前15列\n", ncol(chip_data$soft_platform)))
    }

    datatable(
      preview_df,
      options = list(
        dom = 't',
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "400px",
        columnDefs = list(list(
          className = 'dt-center',
          targets = "_all"
        ))
      ),
      rownames = FALSE,
      filter = 'none',
      escape = FALSE  # 允许HTML渲染（用于显示省略号）
    ) %>%
      formatStyle(columns = 1:ncol(preview_df),
                  fontSize = '85%',
                  maxWidth = '200px',
                  overflow = 'hidden',
                  textOverflow = 'ellipsis')
  })

  # 探针-基因映射预览
  output$chip_probe_mapping_preview <- renderDT({
    req(chip_data$probe_mapping)

    # 显示前10行
    preview_df <- head(chip_data$probe_mapping, 10)

    datatable(
      preview_df,
      options = list(dom = 't', paging = FALSE),
      rownames = FALSE,
      colnames = c("探针ID", "基因符号")
    ) %>%
      formatStyle(columns = c("探针ID", "基因符号"), fontSize = '90%')
  })

  # ============================================
  # 🆕 探针注释与合并工作流 - 服务端逻辑
  # ============================================

  # 步骤2: 渲染SOFT文件ID列选择器
  output$chip_soft_id_column_ui <- renderUI({
    req(chip_data$soft_platform)

    soft_cols <- colnames(chip_data$soft_platform)

    tagList(
      selectInput("chip_soft_id_col",
                  "选择ID列（必须与Series Matrix的探针ID匹配）",
                  choices = c("", soft_cols),  # 添加空选项
                  selected = ""),
      helpText("💡 提示：根据上方SOFT文件预览，选择包含探针ID的列（通常为'ID'列）", style = "color: #ff9800;")
    )
  })

  # 步骤3: 渲染SOFT文件基因列选择器
  output$chip_soft_gene_column_ui <- renderUI({
    req(chip_data$soft_platform)

    soft_cols <- colnames(chip_data$soft_platform)

    tagList(
      selectInput("chip_soft_gene_col",
                  "选择基因列（包含基因符号的列）",
                  choices = c("", soft_cols),  # ✅ 添加空选项
                  selected = ""),  # ✅ 默认不选中，强制用户手动选择
      helpText("⚠️ 请根据上方SOFT文件预览和右侧数据示例，手动选择包含基因符号的列（如 GENE_SYMBOL）。",
               style = "color: #ff9800; font-weight: bold;")
    )
  })

  # 步骤4: 显示选中基因列的实际示例数据
  output$chip_gene_column_examples_ui <- renderUI({
    req(chip_data$soft_platform)
    req(input$chip_soft_gene_col)

    gene_col <- input$chip_soft_gene_col
    examples <- head(chip_data$soft_platform[[gene_col]], 5)

    tagList(
      h6("📋 实际数据示例（前5行）:", style = "color: #667eea;"),
      div(
        style = "background: #f8f9fa; padding: 10px; border-radius: 5px; font-family: monospace; font-size: 11px;",
        for (i in seq_along(examples)) {
          tags$div(
            sprintf("%d. %s", i, as.character(examples[i]))
          )
        }
      )
    )
  })

  # 步骤4: 正则表达式测试
  output$chip_regex_test_result_ui <- renderUI({
    req(input$chip_test_regex)

    # 获取测试数据
    example_text <- if (!is.null(input$chip_gene_extract_example) && input$chip_gene_extract_example != "") {
      input$chip_gene_extract_example
    } else {
      # 使用实际数据示例
      req(input$chip_soft_gene_col)
      examples <- head(chip_data$soft_platform[[input$chip_soft_gene_col]], 3)
      paste(examples, collapse = "\n")
    }

    regex_pattern <- input$chip_gene_regex %||% "[A-Z][A-Z0-9]+"

    # 测试提取
    tryCatch({
      matches <- gregexpr(regex_pattern, example_text, perl = TRUE)
      extracted <- regmatches(example_text, matches)

      # 提取所有匹配项
      all_matches <- unique(unlist(extracted))
      all_matches <- all_matches[all_matches != ""]

      if (length(all_matches) == 0) {
        div(
          class = "alert alert-warning",
          h6("⚠️ 未找到匹配"),
          p("当前正则表达式无法从示例文本中提取基因符号。"),
          p(sprintf("正则表达式: %s", regex_pattern))
        )
      } else {
        tagList(
          div(
            class = "alert alert-success",
            h6("✅ 提取成功", style = "color: #28a745;"),
            p(sprintf("找到 %d 个匹配:", length(all_matches))),
            div(
              style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
              tags$ul(style = "margin: 0; padding-left: 20px;",
                tagList(lapply(all_matches, function(m) {
                  tags$li(style = "margin: 5px 0; font-family: monospace; color: #667eea;",
                          sprintf("<strong>%s</strong>", m))
                }))
              )
            )
          ),
          div(
            style = "margin-top: 10px;",
            h6("📊 提取统计:", style = "color: #667eea;"),
            tags$table(
              class = "table table-striped",
              tags$tbody(
                tags$tr(
                  tags$td("原始文本长度"),
                  tags$td(nchar(example_text))
                ),
                tags$tr(
                  tags$td("提取数量"),
                  tags$td(length(all_matches))
                ),
                tags$tr(
                  tags$td("平均长度"),
                  tags$td(round(mean(nchar(all_matches)), 1))
                )
              )
            )
          )
        )
      }
    }, error = function(e) {
      div(
        class = "alert alert-danger",
        h6("❌ 正则表达式错误"),
        p(e$message)
      )
    })
  })

  # 正则预设按钮逻辑
  observeEvent(input$chip_regex_preset1, {
    updateTextInput(session, "chip_gene_regex", value = "[A-Z]+")
  })

  observeEvent(input$chip_regex_preset2, {
    updateTextInput(session, "chip_gene_regex", value = "[A-Z][A-Z0-9]+")
  })

  observeEvent(input$chip_regex_preset3, {
    updateTextInput(session, "chip_gene_regex", value("\\(([^)]+)\\)"))
  })

  # 步骤5: 预览合并结果
  observeEvent(input$chip_preview_merge, {
    req(chip_data$series_matrix)
    req(chip_data$soft_platform)
    req(input$chip_soft_id_col)
    req(input$chip_soft_gene_col)

    # 转换Series Matrix行名为列
    series_df <- as.data.frame(chip_data$series_matrix)
    if (input$chip_convert_rownames %||% TRUE) {
      series_df <- data.frame(ProbeID = rownames(series_df), series_df, row.names = NULL)
    }

    # 准备SOFT数据
    soft_df <- chip_data$soft_platform[, c(input$chip_soft_id_col, input$chip_soft_gene_col)]
    colnames(soft_df) <- c("ID", "GeneSymbol")

    # 合并数据 - 简化版预览
    merged_df <- merge(series_df, soft_df, by.x = "ProbeID", by.y = "ID", all.x = TRUE)

    # 重新排序列：ProbeID | GeneSymbol | 样本列
    sample_cols <- colnames(merged_df)[!colnames(merged_df) %in% c("ProbeID", "GeneSymbol", "ID")]
    merged_df <- merged_df[, c("ProbeID", "GeneSymbol", sample_cols)]

    # 保存预览结果
    chip_data$merged_preview <- head(merged_df, 5)

    showNotification("✅ 合并预览已生成", type = "message")
  })

  # 显示合并预览表格
  output$chip_merge_preview_table <- renderDT({
    req(chip_data$merged_preview)

    datatable(
      chip_data$merged_preview,
      options = list(
        dom = 't',
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames = FALSE,
      filter = 'none'
    ) %>%
      formatStyle(columns = c("ProbeID", "Gene"), fontSize = '90%')
  })

  # 步骤5: 应用合并配置
  observeEvent(input$chip_apply_merge, {
    req(chip_data$series_matrix)
    req(chip_data$soft_platform)
    req(input$chip_soft_id_col)
    req(input$chip_soft_gene_col)

    showNotification("🔄 正在应用配置并生成最终矩阵...", type = "message")

    # 转换Series Matrix行名为列
    series_df <- as.data.frame(chip_data$series_matrix)
    if (input$chip_convert_rownames %||% TRUE) {
      series_df <- data.frame(ProbeID = rownames(series_df), series_df, row.names = NULL)
    }

    cat(sprintf("🔍 Series探针ID示例: %s\n", paste(head(series_df$ProbeID, 3), collapse = ", ")))

    # 🔍 智能检测SOFT文件中的探针ID列
    user_id_col <- input$chip_soft_id_col
    user_id_sample <- head(chip_data$soft_platform[[user_id_col]][!is.na(chip_data$soft_platform[[user_id_col]])], 10)

    cat(sprintf("🔍 用户选择的ID列 '%s' 示例: %s\n", user_id_col, paste(head(user_id_sample, 3), collapse = ", ")))

    # 检查用户选择的ID列是否真的包含探针ID
    # 支持多种探针ID格式：
    # 1. Affymetrix: 12345_at, 1234567_x_at, 1234567_at
    # 2. Agilent: A_23_P100001, CN_123456
    # 3. Illumina: ILMN_123456
    is_probe_format <- any(grepl(".*_.*_at$|.*_at$|at$", user_id_sample)) ||  # Affymetrix
                       any(grepl("^[A-Z]+_\\d+_", user_id_sample)) ||           # Agilent/Illumina前缀
                       any(grepl("^CN_", user_id_sample))                         # 其他常见格式

    # 准备SOFT数据
    if (is_probe_format) {
      cat("✅ 用户选择的ID列包含探针格式\n")
      soft_df <- chip_data$soft_platform[, c(user_id_col, input$chip_soft_gene_col)]
      colnames(soft_df) <- c("ID", "Gene_Raw")
    } else {
      cat("⚠️ 用户选择的ID列不包含探针格式，尝试使用SOFT行名\n")

      # 尝试使用SOFT文件的行名作为探针ID
      if (!is.null(rownames(chip_data$soft_platform))) {
        # 检查行名是否匹配探针ID
        rowname_sample <- head(rownames(chip_data$soft_platform), 10)
        cat(sprintf("🔍 SOFT行名示例: %s\n", paste(rowname_sample, collapse = ", ")))

        # 创建临时数据框使用行名
        soft_df <- data.frame(
          ID = rownames(chip_data$soft_platform),
          Gene_Raw = chip_data$soft_platform[[input$chip_soft_gene_col]],
          stringsAsFactors = FALSE
        )
      } else {
        cat("❌ SOFT文件没有行名，无法自动检测\n")
        showNotification("⚠️ 请选择包含探针ID的列（如ID、SPOT_ID等）", type = "warning", duration = 10)
        return()
      }
    }

    # 🆕 智能判断是否需要正则提取
    # 检查基因列的内容，如果已经是纯符号格式，不需要提取
    gene_sample <- head(soft_df$Gene_Raw[!is.na(soft_df$Gene_Raw) & soft_df$Gene_Raw != ""], 10)

    cat(sprintf("🔍 用户选择的基因列 '%s' 示例: %s\n",
                input$chip_soft_gene_col, paste(head(gene_sample, 3), collapse = ", ")))

    # 判断是否是纯数字ID
    is_numeric_id <- all(grepl("^[0-9]+$", gene_sample))

    if (is_numeric_id) {
      cat("⚠️ 用户选择的基因列包含数字ID而非基因符号！\n")
      cat("💡 建议：请选择包含基因符号的列（如GENE_SYMBOL）\n")

      # 尝试自动查找真正的基因符号列
      soft_cols <- colnames(chip_data$soft_platform)
      cat(sprintf("🔍 SOFT文件所有列名: %s\n", paste(soft_cols, collapse = ", ")))

      possible_symbol_cols <- c("GENE_SYMBOL", "SYMBOL", "GENE_NAME", "NAME", "DESCRIPTION")

      for (col in possible_symbol_cols) {
        if (col %in% soft_cols && col != input$chip_soft_gene_col) {
          test_data <- head(chip_data$soft_platform[[col]][!is.na(chip_data$soft_platform[[col]])], 10)
          cat(sprintf("🔍 检查列 '%s': 示例=%s\n", col, paste(head(test_data, 3), collapse = ", ")))

          if (!all(grepl("^[0-9]+$", test_data))) {
            cat(sprintf("✅ 自动检测到基因符号列: %s\n", col))
            cat(sprintf("   示例: %s\n", paste(head(test_data, 3), collapse = ", ")))

            # 🔧 修复：创建新的soft_df，直接包含正确的列
            soft_df <- chip_data$soft_platform[, c(user_id_col, col)]
            colnames(soft_df) <- c("ID", "Gene_Raw")
            cat("✅ 已重新创建soft_df，使用正确的基因符号列\n")
            cat(sprintf("✅ 验证：soft_df的Gene_Raw列示例: %s\n", paste(head(soft_df$Gene_Raw[!is.na(soft_df$Gene_Raw)], 3), collapse = ", ")))
            break
          }
        }
      }
    }

    # 重新检查Gene_Raw
    gene_sample <- head(soft_df$Gene_Raw[!is.na(soft_df$Gene_Raw) & soft_df$Gene_Raw != ""], 10)

    # 判断是否已经是纯基因符号（不包含额外文本）
    # 典型基因符号：TP53, EGFR, BRCA1, MYC (大写字母+数字，无空格，无逗号)
    is_pure_symbol <- all(grepl("^[A-Z][A-Z0-9]{1,15}$", gene_sample))

    if (is_pure_symbol) {
      cat("✅ 基因列已是纯符号格式，直接使用，无需正则提取\n")
      cat(sprintf("   示例: %s\n", paste(head(gene_sample, 3), collapse = ", ")))
      # 直接使用原列作为基因符号
      soft_df$GeneSymbol <- soft_df$Gene_Raw
      cat(sprintf("✅ GeneSymbol列已创建，示例: %s\n", paste(head(soft_df$GeneSymbol[!is.na(soft_df$GeneSymbol)], 3), collapse = ", ")))
    } else {
      cat("📋 基因列包含额外文本，应用正则提取\n")
      cat(sprintf("   原始示例: %s\n", paste(head(gene_sample, 2), collapse = ", ")))

      # 应用正则提取
      regex_pattern <- input$chip_gene_regex %||% "[A-Z][A-Z0-9]+"

      soft_df$GeneSymbol <- sapply(soft_df$Gene_Raw, function(x) {
        matches <- regmatches(x, gregexpr(regex_pattern, as.character(x), perl = TRUE))
        if (length(matches[[1]]) > 0) {
          matches[[1]][1]  # 取第一个匹配
        } else {
          NA
        }
      })

      extracted_sample <- head(soft_df$GeneSymbol[!is.na(soft_df$GeneSymbol)], 5)
      cat(sprintf("   提取示例: %s\n", paste(extracted_sample, collapse = ", ")))
    }

    # 🔧 尝试获取Entrez Gene ID（支持用户指定或自动检测）
    # 🔧 通用函数：清理EntrezID（移除非数字字符）
    clean_entrez_id <- function(entrez_str) {
      if (is.na(entrez_str) || is.null(entrez_str) || entrez_str == "") {
        return(NA)
      }
      # 转换为字符
      entrez_str <- as.character(entrez_str)
      # 移除所有非数字字符（保留数字0-9）
      cleaned <- gsub("[^0-9]", "", entrez_str)
      # 如果清理后为空，返回NA
      if (cleaned == "" || is.na(cleaned)) {
        return(NA)
      }
      return(cleaned)
    }

    # 方法1: 用户手动指定EntrezID列
    if (!is.null(input$chip_soft_entrez_col) && input$chip_soft_entrez_col != "" && input$chip_soft_entrez_col != "自动检测") {
      user_entrez_col <- input$chip_soft_entrez_col
      if (user_entrez_col %in% colnames(chip_data$soft_platform)) {
        raw_entrez_ids <- as.character(chip_data$soft_platform[[user_entrez_col]][match(soft_df$ID, chip_data$soft_platform[[input$chip_soft_id_col]])])

        # 🔧 清理EntrezID（移除非数字字符）
        entrez_gene_ids <- sapply(raw_entrez_ids, clean_entrez_id)

        soft_df$EntrezID <- entrez_gene_ids

        # 统计清理情况
        na_count <- sum(is.na(entrez_gene_ids))
        valid_count <- sum(!is.na(entrez_gene_ids))
        entrez_sample <- head(entrez_gene_ids[!is.na(entrez_gene_ids)], 3)

        cat(sprintf("✅ 用户指定的EntrezID列 '%s' 已添加\n", user_entrez_col))
        cat(sprintf("   有效ID数: %d, NA数: %d\n", valid_count, na_count))
        cat(sprintf("   示例: %s\n", paste(entrez_sample, collapse = ", ")))

        if (na_count > 0) {
          cat(sprintf("   ⚠️ %d个ID被清理为NA（包含非数字字符）\n", na_count))
        }
      }
    }
    # 方法2: 自动检测EntrezID列（如果用户选择"自动检测"）
    else if (!is.null(input$chip_soft_entrez_col) && input$chip_soft_entrez_col == "自动检测") {
      # 常见的EntrezID列名（按优先级排序）
      possible_entrez_cols <- c("ENTREZ_GENE_ID", "ENTREZID", "EntrezID", "GeneID",
                                "GENE_ID", "ENTREZ_GENE", "GENE", "Entrez Gene ID",
                                "Entrez", "ENTREZ")

      found_entrez_col <- NULL
      for (col in possible_entrez_cols) {
        if (col %in% colnames(chip_data$soft_platform)) {
          # 检查该列是否包含数字ID（清理后）
          test_values <- head(chip_data$soft_platform[[col]], 100)
          test_values <- test_values[!is.na(test_values) & test_values != ""]

          # 🔧 清理测试值（移除非数字字符）
          cleaned_test <- sapply(test_values, clean_entrez_id)
          cleaned_test <- cleaned_test[!is.na(cleaned_test)]

          # 检查是否大部分可以清理为纯数字（EntrezID特征）
          if (length(cleaned_test) > 0) {
            is_numeric_id <- sum(grepl("^[0-9]+$", cleaned_test)) / length(cleaned_test) > 0.8

            if (is_numeric_id) {
              found_entrez_col <- col
              cat(sprintf("✅ 自动检测到EntrezID列: '%s'\n", col))
              break
            }
          }
        }
      }

      if (!is.null(found_entrez_col)) {
        raw_entrez_ids <- as.character(chip_data$soft_platform[[found_entrez_col]][match(soft_df$ID, chip_data$soft_platform[[input$chip_soft_id_col]])])

        # 🔧 清理EntrezID（移除非数字字符）
        entrez_gene_ids <- sapply(raw_entrez_ids, clean_entrez_id)

        soft_df$EntrezID <- entrez_gene_ids

        # 统计清理情况
        na_count <- sum(is.na(entrez_gene_ids))
        valid_count <- sum(!is.na(entrez_gene_ids))
        entrez_sample <- head(entrez_gene_ids[!is.na(entrez_gene_ids)], 3)

        cat(sprintf("✅ 自动检测的EntrezID列已添加\n", found_entrez_col))
        cat(sprintf("   有效ID数: %d, NA数: %d\n", valid_count, na_count))
        cat(sprintf("   示例: %s\n", paste(entrez_sample, collapse = ", ")))

        if (na_count > 0) {
          cat(sprintf("   ⚠️ %d个ID被清理为NA（包含非数字字符）\n", na_count))
        }
      } else {
        cat("⚠️ 自动检测未找到EntrezID列\n")
      }
    }
    # 方法3: 旧逻辑（仅检查GENE列）- 也添加清理
    else if ("GENE" %in% colnames(chip_data$soft_platform)) {
      raw_entrez_ids <- as.character(chip_data$soft_platform$GENE[match(soft_df$ID, chip_data$soft_platform[[input$chip_soft_id_col]])])

      # 🔧 清理EntrezID
      entrez_gene_ids <- sapply(raw_entrez_ids, clean_entrez_id)

      soft_df$EntrezID <- entrez_gene_ids

      # 统计清理情况
      na_count <- sum(is.na(entrez_gene_ids))
      valid_count <- sum(!is.na(entrez_gene_ids))
      entrez_sample <- head(entrez_gene_ids[!is.na(entrez_gene_ids)], 3)

      cat(sprintf("✅ Entrez Gene ID列（GENE）已添加\n"))
      cat(sprintf("   有效ID数: %d, NA数: %d\n", valid_count, na_count))
      cat(sprintf("   示例: %s\n", paste(entrez_sample, collapse = ", ")))

      if (na_count > 0) {
        cat(sprintf("   ⚠️ %d个ID被清理为NA（包含非数字字符）\n", na_count))
      }
    }

    # 🔧 不再创建Gene列！直接使用GeneSymbol用于去重，EntrezID用于差异分析
    cat("✅ 使用GeneSymbol列进行探针去重\n")
    cat("✅ 使用EntrezID列作为差异分析结果的ID\n")

    # 合并数据 - 保留所有列
    merged_df <- merge(series_df, soft_df, by.x = "ProbeID", by.y = "ID", all.x = TRUE)

    # 🔧 智能识别样本列（必须是数值型，排除ID列和基因列）
    # 规则：样本列 = 数值型列 + 不在排除列表中
    exclude_cols <- c("ProbeID", "GeneSymbol", "EntrezID", "Gene_Raw", "ID")

    # 识别数值型列（样本列）
    sample_cols <- character(0)
    for (col in colnames(merged_df)) {
      if (!(col %in% exclude_cols)) {
        # 检查是否为数值型
        if (is.numeric(merged_df[[col]])) {
          sample_cols <- c(sample_cols, col)
        }
      }
    }

    cat(sprintf("🔍 识别到 %d 个样本列: %s\n", length(sample_cols),
                paste(head(sample_cols, 5), collapse = ", ")))

    # 重新排序列：ProbeID | GeneSymbol | EntrezID | 样本列
    if ("EntrezID" %in% colnames(merged_df)) {
      # 重新排列列顺序
      merged_df <- merged_df[, c("ProbeID", "GeneSymbol", "EntrezID", sample_cols)]

      # 🔧 强制确保EntrezID是字符型
      merged_df$EntrezID <- as.character(merged_df$EntrezID)
      cat("✅ EntrezID列已强制转换为字符型\n")
    } else {
      # 如果没有EntrezID，只保留GeneSymbol
      merged_df <- merged_df[, c("ProbeID", "GeneSymbol", sample_cols)]
      cat("⚠️ 未找到EntrezID列\n")
    }

    # 🔍 诊断：显示合并后的列结构
    cat(sprintf("🔍 合并后矩阵结构: %d 行 × %d 列\n", nrow(merged_df), ncol(merged_df)))
    cat(sprintf("🔍 列名: %s\n", paste(colnames(merged_df), collapse = ", ")))
    cat(sprintf("🔍 前3列: %s\n", paste(head(colnames(merged_df), 3), collapse = ", ")))
    cat(sprintf("🔍 后3列: %s\n", paste(tail(colnames(merged_df), 3), collapse = ", ")))

    # 显示各列的类型
    cat("🔍 列类型:\n")
    for (i in 1:min(5, ncol(merged_df))) {
      cat(sprintf("   %s: %s\n", colnames(merged_df)[i], class(merged_df[[i]])))
    }

    # 确认EntrezID列的位置和类型
    if ("EntrezID" %in% colnames(merged_df)) {
      entrez_col_idx <- which(colnames(merged_df) == "EntrezID")
      cat(sprintf("🔍 EntrezID列位置: 第%d列，类型: %s\n", entrez_col_idx, class(merged_df$EntrezID)))
      cat(sprintf("🔍 EntrezID示例: %s\n", paste(head(merged_df$EntrezID[!is.na(merged_df$EntrezID)], 3), collapse = ", ")))
    }

    # 保存最终合并的矩阵
    chip_data$merged_matrix <- merged_df

    # 🆕 检查数据列内容
    gene_symbol_sample <- head(merged_df$GeneSymbol[!is.na(merged_df$GeneSymbol)], 5)
    cat(sprintf("✅ GeneSymbol列包含基因符号（示例: %s）\n",
                paste(gene_symbol_sample, 3), collapse = ", "))

    if ("EntrezID" %in% colnames(merged_df)) {
      entrez_sample <- head(merged_df$EntrezID[!is.na(merged_df$EntrezID)], 5)
      cat(sprintf("✅ EntrezID列包含Entrez Gene ID（示例: %s）\n",
                  paste(entrez_sample, collapse = ", ")))
    }

    # 统计信息
    n_total <- nrow(merged_df)
    n_annotated <- sum(!is.na(merged_df$GeneSymbol))
    annotation_rate <- n_annotated / n_total * 100

    # 🔍 诊断信息：检查GeneSymbol列匹配情况
    na_count <- sum(is.na(merged_df$GeneSymbol) | merged_df$GeneSymbol == "")
    cat(sprintf("🔍 诊断: %d个探针的GeneSymbol为NA或空 (%.1f%%)\n", na_count, na_count/n_total*100))

    # 检查ID匹配情况
    series_probes <- series_df$ProbeID
    soft_ids <- soft_df$ID
    matched_probes <- sum(series_probes %in% soft_ids)
    cat(sprintf("🔍 ID匹配: %d / %d 个Series探针在SOFT文件中找到 (%.1f%%)\n",
                matched_probes, length(series_probes), matched_probes/length(series_probes)*100))

    # 如果匹配率太低，显示警告
    if (matched_probes / length(series_probes) < 0.5) {
      cat("⚠️ 警告：ID匹配率低于50%，请检查ID列选择是否正确！\n")
      cat(sprintf("   Series探针示例: %s\n", paste(head(series_probes, 3), collapse = ", ")))
      cat(sprintf("   SOFT ID示例: %s\n", paste(head(soft_ids, 3), collapse = ", ")))
    }

    cat(sprintf("✅ 合并完成: %d 个探针, %d 个已注释 (%.1f%%)\n",
                n_total, n_annotated, annotation_rate))

    showNotification(
      sprintf("✅ 合并完成！%d / %d 个探针已注释 (%.1f%%)",
              n_annotated, n_total, annotation_rate),
      type = "message",
      duration = 10
    )
  })

  # 显示最终合并矩阵的预览（在应用合并后）
  output$chip_final_matrix_ui <- renderUI({
    req(chip_data$merged_matrix)

    tagList(
      h5("📊 最终表达矩阵（前5行）", style = "color: #28a745;"),
      DTOutput("chip_final_matrix_table")
    )
  })

  output$chip_final_matrix_table <- renderDT({
    req(chip_data$merged_matrix)

    # 显示前5行和前10列（避免表格过宽）
    preview_df <- head(chip_data$merged_matrix, 5)
    if (ncol(preview_df) > 12) {
      preview_df <- preview_df[, 1:12]  # ProbeID + Gene + 前10个样本
    }

    # 检查是否有ProbeID和Gene列
    has_probe_id <- "ProbeID" %in% colnames(preview_df)
    has_gene <- "Gene" %in% colnames(preview_df)

    # 创建datatable
    dt <- datatable(
      preview_df,
      options = list(
        dom = 't',
        paging = FALSE,
        scrollX = TRUE,
        columnDefs = list(list(
          className = 'dt-center',
          targets = "_all"
        ))
      ),
      rownames = FALSE,
      filter = 'none'
    )

    # 只对存在的列应用样式
    if (has_probe_id && has_gene) {
      dt <- dt %>%
        formatStyle(columns = c("ProbeID", "Gene"),
                    backgroundColor = '#e8f4f8',
                    fontWeight = 'bold')
    }

    dt
  })

  # ============================================
  # 🆕 数据预处理与探针去重 - Server逻辑
  # ============================================

  # 去重前统计显示
  output$chip_before_dedupe_stats <- renderUI({
    req(chip_data$merged_matrix)

    # 统计探针数量
    n_probes <- nrow(chip_data$merged_matrix)

    # 识别数值列（样本数据）
    numeric_cols <- sapply(chip_data$merged_matrix, function(x) is.numeric(x))
    n_samples <- sum(numeric_cols)

    # 统计探针-基因映射情况
    n_with_gene <- sum(!is.na(chip_data$merged_matrix$GeneSymbol))

    # 计算一个基因对应多个探针的情况
    gene_counts <- table(chip_data$merged_matrix$GeneSymbol)
    gene_counts <- gene_counts[names(gene_counts) != ""]  # 移除空基因名
    n_multi_probes <- sum(gene_counts > 1)

    div(
      style = "font-size: 12px;",
      tags$ul(style = "padding-left: 15px; margin: 5px 0;",
        tags$li(sprintf("总探针数: %d", n_probes)),
        tags$li(sprintf("样本数: %d", n_samples)),
        tags$li(sprintf("有基因注释的探针: %d (%.1f%%)", n_with_gene, n_with_gene/n_probes*100)),
        tags$li(sprintf("一因多探针的基因数: %d", n_multi_probes))
      )
    )
  })

  # 步骤1: 数据预处理
  observeEvent(input$chip_preprocess_data, {
    req(chip_data$merged_matrix)

    showNotification("🔄 正在进行数据预处理...", type = "message")

    tryCatch({
      # 🔧 修复：智能提取表达数据（只保留数值列）
      merged_df <- chip_data$merged_matrix

      # 识别数值列
      numeric_cols <- sapply(merged_df, function(x) is.numeric(x))

      # 排除ProbeID和Gene列（它们在前两列）
      # 但要确保只使用数值列进行计算
      expr_cols <- which(numeric_cols)

      # 提取表达数据
      expr_matrix <- as.matrix(merged_df[, expr_cols, drop = FALSE])

      # 🔧 修复：使用ProbeID作为行名，而不是Gene（这样才能在去重时正确匹配）
      if ("ProbeID" %in% colnames(merged_df)) {
        rownames(expr_matrix) <- merged_df$ProbeID
        cat("✅ 表达矩阵行名使用ProbeID\n")
      } else if ("Gene" %in% colnames(merged_df)) {
        rownames(expr_matrix) <- merged_df$Gene
        cat("⚠️ 警告：使用Gene作为行名（ProbeID列不存在）\n")
      } else {
        rownames(expr_matrix) <- rownames(merged_df)
      }

      cat(sprintf("✅ 提取表达数据: %d 探针 × %d 样本\n",
                  nrow(expr_matrix), ncol(expr_matrix)))

      # 保存原始数据用于对比
      chip_data$expr_before_preprocess <- expr_matrix

      # 1. 自动log2转换判断
      if (input$chip_auto_log2 %||% TRUE) {
        ex <- expr_matrix
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
        LogC <- (qx[5] > 100) ||
                (qx[6] - qx[1] > 50 && qx[2] > 0) ||
                (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

        if (LogC) {
          ex[which(ex <= 0)] <- NaN
          expr_matrix <- log2(ex)
          chip_data$log2_performed <- TRUE
          cat("✅ log2转换已完成\n")
        } else {
          chip_data$log2_performed <- FALSE
          cat("ℹ️ 不需要log2转换\n")
        }
      }

      # 2. limma标准化
      if (input$chip_normalize_data %||% TRUE) {
        library(limma)
        expr_matrix <- normalizeBetweenArrays(expr_matrix)
        chip_data$normalize_performed <- TRUE
        cat("✅ limma标准化已完成\n")
      }

      # 保存处理后的数据
      chip_data$expr_preprocessed <- expr_matrix

      # 生成结果报告
      chip_data$preprocess_report <- list(
        log2_performed = chip_data$log2_performed,
        normalize_performed = chip_data$normalize_performed,
        n_probes = nrow(expr_matrix),
        n_samples = ncol(expr_matrix),
        data_range = range(expr_matrix, na.rm = TRUE)
      )

      # 🆕 生成箱线图对比（矫正前后）
      tryCatch({
        library(ggplot2)
        library(reshape2)

        # 矫正前的数据 - 计算每个样本的统计值用于箱线图
        expr_before <- chip_data$expr_before_preprocess
        # 转置矩阵：行为样本，列为探针
        df_before <- as.data.frame(t(expr_before))
        df_before$Sample <- rownames(df_before)
        # 将数据从宽格式转换为长格式（使用reshape2的melt）
        df_before_long <- melt(df_before, id.vars = "Sample", variable.name = "Probe", value.name = "Expression")
        df_before_long$Stage <- "Before"

        # 矫正后的数据
        df_after <- as.data.frame(t(expr_matrix))
        df_after$Sample <- rownames(df_after)
        df_after_long <- melt(df_after, id.vars = "Sample", variable.name = "Probe", value.name = "Expression")
        df_after_long$Stage <- "After"

        # 合并数据
        df_combined <- rbind(df_before_long, df_after_long)

        # 绘制箱线图
        p <- ggplot(df_combined, aes(x = Sample, y = Expression, fill = Stage)) +
          geom_boxplot() +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            legend.position = "top",
            plot.title = element_text(hjust = 0.5)
          ) +
          labs(
            title = "数据预处理前后对比",
            x = "样本",
            y = "表达值",
            fill = "阶段"
          ) +
          scale_fill_manual(values = c("Before" = "#E69F00", "After" = "#009E73"))

        chip_data$preprocess_boxplot <- p
        cat("✅ 箱线图已生成\n")
      }, error = function(e) {
        cat(sprintf("⚠️ 箱线图生成失败: %s\n", e$message))
      })

      showNotification("✅ 数据预处理完成！", type = "message", duration = 5)

    }, error = function(e) {
      showNotification(sprintf("❌ 预处理失败: %s", e$message), type = "error")
    })
  })

  # 预处理结果显示
  output$chip_preprocess_result_ui <- renderUI({
    req(chip_data$preprocess_report)
    req(chip_data$expr_preprocessed)

    report <- chip_data$preprocess_report
    expr_matrix <- chip_data$expr_preprocessed

    # 计算额外统计信息
    expr_mean <- mean(expr_matrix, na.rm = TRUE)
    expr_median <- median(expr_matrix, na.rm = TRUE)
    expr_sd <- sd(expr_matrix, na.rm = TRUE)

    tagList(
      h5("✅ 预处理完成", style = "color: #28a745; font-weight: bold;"),
      br(),

      # 基本统计信息
      wellPanel(
        style = "background: #f8f9fa; border: 1px solid #dee2e6;",
        h6("📊 基本统计", style = "color: #495057; margin-top: 0;"),
        tags$table(
          class = "table table-sm table-striped",
          style = "margin-bottom: 0;",
          tags$thead(
            tags$tr(
              tags$th("项目", style = "width: 50%;"),
              tags$th("值")
            )
          ),
          tags$tbody(
            tags$tr(
              tags$td("log2转换"),
              tags$td({
                if (report$log2_performed)
                  tags$span("✅ 是", class = "badge badge-success")
                else
                  tags$span("❌ 否", class = "badge badge-secondary")
              })
            ),
            tags$tr(
              tags$td("limma标准化 (quantile)"),
              tags$td({
                if (report$normalize_performed)
                  tags$span("✅ 是", class = "badge badge-success")
                else
                  tags$span("❌ 否", class = "badge badge-secondary")
              })
            ),
            tags$tr(
              tags$td(tags$strong("探针数")),
              tags$td(sprintf("%d", report$n_probes))
            ),
            tags$tr(
              tags$td(tags$strong("样本数")),
              tags$td(sprintf("%d", report$n_samples))
            ),
            tags$tr(
              tags$td("数据范围"),
              tags$td(sprintf("%.3f ~ %.3f", report$data_range[1], report$data_range[2]))
            ),
            tags$tr(
              tags$td("平均值"),
              tags$td(sprintf("%.3f", expr_mean))
            ),
            tags$tr(
              tags$td("中位数"),
              tags$td(sprintf("%.3f", expr_median))
            ),
            tags$tr(
              tags$td("标准差"),
              tags$td(sprintf("%.3f", expr_sd))
            )
          )
        )
      ),
      br(),

      # 箱线图对比
      wellPanel(
        style = "background: white; border: 1px solid #dee2e6;",
        h6("📊 箱线图对比（矫正前后）", style = "color: #ff9800; margin-top: 0;"),
        plotOutput("chip_preprocess_boxplot", height = "500px")
      )
    )
  })

  # 🆕 渲染箱线图
  output$chip_preprocess_boxplot <- renderPlot({
    req(chip_data$preprocess_boxplot)
    chip_data$preprocess_boxplot
  })

  # ============================================
  # 🆕 批次矫正 - Server逻辑
  # ============================================

  # 自动检测批次
  observeEvent(input$chip_detect_batch, {
    req(chip_data$expr_preprocessed)
    req(input$chip_batch_pattern)

    pattern <- input$chip_batch_pattern
    sample_names <- colnames(chip_data$expr_preprocessed)

    # 从样本名中提取批次信息
    batches <- sapply(sample_names, function(x) {
      match <- regmatches(x, regexpr(paste0(pattern, "\\w+"), x, ignore.case = TRUE))
      if (length(match) > 0) {
        return(match)
      } else {
        return("Unknown")
      }
    })

    # 保存批次信息
    chip_data$batch_info <- batches
    chip_data$batch_method <- input$chip_batch_method %||% "limma"

    # 统计批次分布
    batch_table <- table(batches)

    cat(sprintf("✅ 检测到 %d 个批次:\n", length(batch_table)))
    for (i in seq_along(batch_table)) {
      cat(sprintf("   - %s: %d 个样本\n", names(batch_table)[i], batch_table[i]))
    }

    showNotification(
      sprintf("✅ 检测到 %d 个批次", length(batch_table)),
      type = "message"
    )
  })

  # 解析手动指定的批次
  observeEvent(input$chip_batch_manual, {
    req(input$chip_batch_manual)

    tryCatch({
      lines <- strsplit(input$chip_batch_manual, "\n")[[1]]
      lines <- lines[lines != "" & !grepl("^\\s*$", lines)]

      batch_mapping <- list()
      for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) {
          batch_mapping[[trimws(parts[1])]] <- trimws(parts[2])
        }
      }

      chip_data$batch_manual_mapping <- batch_mapping

      sample_names <- colnames(chip_data$expr_preprocessed)
      batches <- sapply(sample_names, function(x) {
        if (x %in% names(batch_mapping)) {
          return(batch_mapping[[x]])
        } else {
          return("Unknown")
        }
      })

      chip_data$batch_info <- batches

      batch_table <- table(batches)

      cat(sprintf("✅ 手动指定批次: %d 个批次\n", length(batch_table)))

      showNotification(
        sprintf("✅ 手动指定批次: %d 个批次", length(batch_table)),
        type = "message"
      )
    }, error = function(e) {
      showNotification(sprintf("❌ 解析批次信息失败: %s", e$message), type = "error")
    })
  })

  # 显示批次信息
  output$chip_batch_info_ui <- renderUI({
    req(chip_data$batch_info)

    batches <- chip_data$batch_info
    batch_table <- table(batches)

    tagList(
      tags$table(
        class = "table table-striped",
        tags$thead(
          tags$tr(
            tags$th("批次"),
            tags$th("样本数")
          )
        ),
        tags$tbody(
          lapply(seq_along(batch_table), function(i) {
            tags$tr(
              tags$td(names(batch_table)[i]),
              tags$td(batch_table[i])
            )
          })
        )
      )
    )
  })

  # 执行批次矫正
  observeEvent(input$chip_apply_batch_correct, {
    req(chip_data$expr_preprocessed)
    req(chip_data$batch_info)

    showNotification("🔧 正在进行批次矫正...", type = "message")

    tryCatch({
      expr_matrix <- chip_data$expr_preprocessed
      batches <- chip_data$batch_info
      method <- input$chip_batch_method %||% "limma"

      cat(sprintf("📊 批次矫正方法: %s\n", method))

      if (method == "limma") {
        # 使用 limma::removeBatchEffect
        library(limma)

        batch_factor <- factor(batches)
        design <- model.matrix(~1, data = data.frame(batch = batch_factor))

        expr_corrected <- removeBatchEffect(expr_matrix, batch = batch_factor)

        cat("✅ limma::removeBatchEffect 批次矫正完成\n")

      } else if (method == "combat") {
        # 使用 sva::ComBat
        library(sva)

        batch_factor <- factor(batches)

        expr_corrected <- ComBat(
          dat = expr_matrix,
          batch = batch_factor,
          mod = NULL,
          par.prior = TRUE,
          prior.plots = FALSE
        )

        cat("✅ sva::ComBat 批次矫正完成\n")
      }

      # 保存矫正后的数据
      chip_data$expr_batch_corrected <- expr_corrected
      chip_data$batch_correct_method <- method

      # 计算矫正前后的差异
      mean_diff <- mean(abs(expr_matrix - expr_corrected), na.rm = TRUE)

      # 生成报告
      chip_data$batch_correct_report <- list(
        method = method,
        n_batches = length(unique(batches)),
        batch_distribution = table(batches),
        mean_change = mean_diff,
        n_probes = nrow(expr_corrected),
        n_samples = ncol(expr_corrected)
      )

      showNotification(
        sprintf("✅ 批次矫正完成！方法: %s", method),
        type = "message",
        duration = 10
      )

    }, error = function(e) {
      showNotification(
        sprintf("❌ 批次矫正失败: %s", e$message),
        type = "error",
        duration = 10
      )
      cat(sprintf("❌ 批次矫正错误: %s\n", e$message))
    })
  })

  # 显示批次矫正结果
  output$chip_batch_correct_result_ui <- renderUI({
    req(chip_data$batch_correct_report)

    report <- chip_data$batch_correct_report

    method_name <- if (report$method == "limma") {
      "limma::removeBatchEffect"
    } else {
      "sva::ComBat"
    }

    tagList(
      h6("✅ 批次矫正完成", style = "color: #E91E63;"),
      tags$table(
        class = "table table-striped",
        tags$thead(
          tags$tr(
            tags$th("项目"),
            tags$th("值")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td("矫正方法"),
            tags$td(method_name)
          ),
          tags$tr(
            tags$td("批次数"),
            tags$td(report$n_batches)
          ),
          tags$tr(
            tags$td("探针/基因数"),
            tags$td(report$n_probes)
          ),
          tags$tr(
            tags$td("样本数"),
            tags$td(report$n_samples)
          ),
          tags$tr(
            tags$td("平均变化幅度"),
            tags$td(sprintf("%.4f", report$mean_change))
          )
        )
      ),
      br(),
      helpText("💡 提示：批次矫正后的数据已保存，将用于后续的探针去重和差异分析。")
    )
  })

  # 步骤2: 探针去重
  observeEvent(input$chip_dedupe_probes, {
    req(chip_data$expr_preprocessed)

    showNotification("✂️ 正在进行探针去重...", type = "message")

    tryCatch({
      library(dplyr)
      library(tibble)

      # ✅ 优先使用批次矫正后的数据
      if (!is.null(chip_data$expr_batch_corrected)) {
        expr_matrix <- chip_data$expr_batch_corrected
        cat("✅ 使用批次矫正后的数据进行探针去重\n")
      } else {
        expr_matrix <- chip_data$expr_preprocessed
        cat("✅ 使用预处理后的数据进行探针去重\n")
      }

      # 转换为数据框并添加探针ID和基因信息
      expr_df <- as.data.frame(expr_matrix)
      expr_df <- expr_df %>%
        rownames_to_column("ProbeID")

      # 🔧 检查merged_matrix中的可用列
      available_cols <- colnames(chip_data$merged_matrix)
      cat(sprintf("📋 merged_matrix可用列: %s\n", paste(available_cols, collapse = ", ")))

      # 动态选择要合并的列（只选择存在的列）
      merge_cols <- c("ProbeID")
      if ("GeneSymbol" %in% available_cols) {
        merge_cols <- c(merge_cols, "GeneSymbol")
      }
      if ("EntrezID" %in% available_cols) {
        merge_cols <- c(merge_cols, "EntrezID")
      }

      cat(sprintf("📋 将合并列: %s\n", paste(merge_cols, collapse = ", ")))

      # 执行合并
      expr_df <- expr_df %>%
        inner_join(chip_data$merged_matrix[, merge_cols, drop = FALSE], by = "ProbeID")

      # 🔧 清理合并后的EntrezID列（移除非数字字符）
      if ("EntrezID" %in% colnames(expr_df)) {
        clean_entrez_id <- function(entrez_str) {
          if (is.na(entrez_str) || is.null(entrez_str) || entrez_str == "") {
            return(NA)
          }
          entrez_str <- as.character(entrez_str)
          # 移除所有非数字字符（///, //, -, 等）
          cleaned <- gsub("[^0-9]", "", entrez_str)
          if (cleaned == "" || is.na(cleaned)) {
            return(NA)
          }
          return(cleaned)
        }

        # 统计清理前的情况
        na_before <- sum(is.na(expr_df$EntrezID))

        # 应用清理
        expr_df$EntrezID <- sapply(expr_df$EntrezID, clean_entrez_id)

        # 统计清理后的情况
        na_after <- sum(is.na(expr_df$EntrezID))

        cat(sprintf("🔧 EntrezID清理: NA前=%d, NA后=%d, 新增NA=%d\n",
                    na_before, na_after, na_after - na_before))

        if (na_after > na_before) {
          cat(sprintf("⚠️ %d个EntrezID包含非数字字符，已清理为NA\n", na_after - na_before))
        }
      }

      # 🔍 诊断：检查合并后的GeneSymbol列
      cat(sprintf("🔍 去重前: %d 行，GeneSymbol NA=%d, 空字符串=%d\n",
                  nrow(expr_df),
                  sum(is.na(expr_df$GeneSymbol)),
                  sum(expr_df$GeneSymbol == "")))

      # 移除没有基因符号的探针
      expr_df <- expr_df[!is.na(expr_df$GeneSymbol) & expr_df$GeneSymbol != "", ]

      cat(sprintf("🔍 移除NA后: %d 行剩余\n", nrow(expr_df)))

      # 如果移除后没有数据，警告并跳过去重
      if (nrow(expr_df) == 0) {
        cat("⚠️ 警告：移除NA基因后没有数据！请检查合并步骤的GeneSymbol列\n")
        showNotification("❌ 去重失败：所有基因都是NA，请检查ID列和基因列是否匹配", type = "error", duration = 10)
        return()
      }

      # 探针去重：保留表达量最高的探针
      # 使用GeneSymbol去重，但保留EntrezID列供差异分析使用
      # 最终expr_deduped包含：行名=GeneSymbol, 列=EntrezID + 样本数据
      if ("EntrezID" %in% colnames(expr_df)) {
        expr_df <- expr_df %>%
          select(-ProbeID) %>%                    # 去掉探针ID列
          select(GeneSymbol, EntrezID, everything()) %>%  # 基因符号和ID放前面
          mutate(rowMean = rowMeans(.[, -(1:2)])) %>% # 计算平均表达量（排除前两列）
          arrange(desc(rowMean)) %>%               # 按表达量降序排列
          distinct(GeneSymbol, .keep_all = TRUE) %>%     # 按基因符号去重
          select(-rowMean) %>%                     # 删除辅助列（保留EntrezID列！）
          column_to_rownames("GeneSymbol")         # GeneSymbol作为行名，EntrezID保留为列

        cat("✅ 去重使用GeneSymbol，行名=GeneSymbol，EntrezID保留为列供差异分析使用\n")
      } else {
        # 没有EntrezID，使用GeneSymbol
        expr_df <- expr_df %>%
          select(-ProbeID) %>%                    # 去掉探针ID列
          select(GeneSymbol, everything()) %>%   # 基因符号放第一
          mutate(rowMean = rowMeans(.[, -1])) %>% # 计算平均表达量
          arrange(desc(rowMean)) %>%               # 按表达量降序排列
          distinct(GeneSymbol, .keep_all = TRUE) %>%     # 去重
          select(-rowMean) %>%                     # 删除辅助列
          column_to_rownames("GeneSymbol")        # 基因符号变回行名

        cat("✅ 去重使用GeneSymbol，最终结果使用GeneSymbol作为行名\n")
      }

      # 保存去重后的数据
      chip_data$expr_deduped <- expr_df

      # 统计信息
      n_before <- nrow(chip_data$expr_preprocessed)
      n_after <- nrow(expr_df)
      reduction_rate <- (n_before - n_after) / n_before * 100

      chip_data$dedupe_report <- list(
        n_probes_before = n_before,
        n_genes_after = n_after,
        n_removed = n_before - n_after,
        reduction_rate = reduction_rate,
        n_samples = ncol(expr_df)
      )

      cat(sprintf("✅ 探针去重完成: %d 探针 → %d 基因 (%.1f%% 减少)\n",
                  n_before, n_after, reduction_rate))

      showNotification(
        sprintf("✅ 去重完成！%d 探针 → %d 基因", n_before, n_after),
        type = "message",
        duration = 5
      )

    }, error = function(e) {
      showNotification(sprintf("❌ 去重失败: %s", e$message), type = "error")
    })
  })

  # 去重后统计显示
  output$chip_after_dedupe_stats <- renderUI({
    req(chip_data$dedupe_report)

    report <- chip_data$dedupe_report

    div(
      style = "font-size: 12px; color: #28a745;",
      tags$ul(style = "padding-left: 15px; margin: 5px 0;",
        tags$li(sprintf("基因数: %d", report$n_genes_after)),
        tags$li(sprintf("样本数: %d", report$n_samples)),
        tags$li(sprintf("减少探针数: %d (%.1f%%)", report$n_removed, report$reduction_rate))
      )
    )
  })

  # 去重结果显示
  output$chip_dedupe_result_ui <- renderUI({
    req(chip_data$dedupe_report)

    report <- chip_data$dedupe_report

    tagList(
      tags$table(
        class = "table table-striped",
        tags$thead(
          tags$tr(
            tags$th("项目"),
            tags$th("值")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td("去重前探针数"),
            tags$td(report$n_probes_before)
          ),
          tags$tr(
            tags$td("去重后基因数"),
            tags$td(report$n_genes_after)
          ),
          tags$tr(
            tags$td("减少的探针数"),
            tags$td(report$n_removed)
          ),
          tags$tr(
            tags$td("减少比例"),
            tags$td(sprintf("%.1f%%", report$reduction_rate))
          ),
          tags$tr(
            tags$td("样本数"),
            tags$td(report$n_samples)
          )
        )
      ),

      br(),

      h6("📊 去重后表达矩阵预览（前5行 × 前10列）：", style = "color: #666;"),

      # 显示去重后的矩阵预览
      DTOutput("chip_deduped_matrix_preview")
    )
  })

  # 去重后矩阵预览表格
  output$chip_deduped_matrix_preview <- renderDT({
    req(chip_data$expr_deduped)

    preview_df <- head(chip_data$expr_deduped, 5)
    if (ncol(preview_df) > 10) {
      preview_df <- preview_df[, 1:10]
    }

    datatable(
      as.data.frame(preview_df),
      options = list(
        dom = 't',
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames = TRUE,
      filter = 'none'
    ) %>%
      formatStyle(columns = 1:ncol(preview_df), fontSize = '85%')
  })

  # 步骤3: 生成标准格式数据
  observeEvent(input$chip_generate_standard_data, {
    req(chip_data$expr_deduped)

    showNotification("🚀 正在生成标准格式数据...", type = "message")

    tryCatch({
      # 保存为标准格式（可直接用于现有分析模块）
      chip_data$standard_expression <- chip_data$expr_deduped

      # 保存样本名称
      chip_data$sample_names <- colnames(chip_data$expr_deduped)

      # 保存基因名称
      chip_data$gene_names <- rownames(chip_data$expr_deduped)

      # 标记数据已准备好
      chip_data$ready_for_analysis <- TRUE

      cat(sprintf("✅ 标准格式数据已生成: %d 基因 × %d 样本\n",
                  nrow(chip_data$standard_expression),
                  ncol(chip_data$standard_expression)))

      showNotification(
        sprintf("✅ 标准格式数据已生成！%d 基因 × %d 样本",
                nrow(chip_data$standard_expression),
                ncol(chip_data$standard_expression)),
        type = "message",
        duration = 10
      )

    }, error = function(e) {
      showNotification(sprintf("❌ 生成失败: %s", e$message), type = "error")
    })
  })

  # 标准数据摘要显示
  output$chip_standard_data_summary <- renderUI({
    req(chip_data$standard_expression)

    tagList(
      h6("📊 数据摘要", style = "color: #155724;"),
      tags$table(
        class = "table table-striped",
        tags$thead(
          tags$tr(
            tags$th("项目"),
            tags$th("值")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td("基因数"),
            tags$td(nrow(chip_data$standard_expression))
          ),
          tags$tr(
            tags$td("样本数"),
            tags$td(ncol(chip_data$standard_expression))
          ),
          tags$tr(
            tags$td("样本名称"),
            tags$td(paste(head(chip_data$sample_names, 5), collapse = ", "))
          )
        )
      ),
      br(),
      h6("✅ 现在可以进行以下分析：", style = "color: #155724;"),
      tags$ul(style = "padding-left: 20px;",
        tags$li("切换到“差异分析”模块进行limma分析"),
        tags$li("使用样本分组功能设置对照组和处理组"),
        tags$li("进行KEGG和GO富集分析"),
        tags$li("生成火山图和其他可视化")
      )
    )
  })

  # ============================================
  # 原有的分组逻辑继续
  # ============================================

  # 动态生成分组 UI
  output$chip_grouping_ui <- renderUI({
    req(chip_data$series_matrix)

    sample_names <- colnames(chip_data$series_matrix)

    # 如果自动检测到分组
    if (!is.null(chip_data$group_info) && !is.null(chip_data$group_info$pattern_name)) {
      tagList(
        div(
          class = "alert alert-info",
          h5("✅ 自动检测到分组模式"),
          p(sprintf("模式: %s", chip_data$group_info$pattern_name)),
          p(sprintf("对照组: %s",
                    paste(chip_data$group_info$ctrl_samples, collapse = ", "))),
          p(sprintf("处理组: %s",
                    paste(chip_data$group_info$trt_samples, collapse = ", "))),
          checkboxInput("chip_use_auto_groups",
                       "使用自动检测的分组",
                       value = TRUE)
        ),

        conditionalPanel(
          condition = "!input.chip_use_auto_groups",
          h5("手动选择分组:"),
          helpText("💡 提示：可以点击输入框后，直接粘贴样本名称（用逗号或空格分隔）"),
          fluidRow(
            column(6,
              selectizeInput("chip_ctrl_samples",
                         "对照组样本:",
                         choices = sample_names,
                         multiple = TRUE,
                         options = list(create = TRUE))
            ),
            column(6,
              selectizeInput("chip_trt_samples",
                         "处理组样本:",
                         choices = sample_names,
                         multiple = TRUE,
                         options = list(create = TRUE))
            )
          )
        )
      )
    } else {
      # 手动分组
      tagList(
        div(
          class = "alert alert-warning",
          h5("⚠️ 未能自动检测分组模式"),
          p("请使用上方的快速粘贴功能，或手动指定对照组和处理组样本。")
        ),
        helpText("💡 提示：可以点击输入框后，直接粘贴样本名称（用逗号或空格分隔）"),
        fluidRow(
          column(6,
            selectizeInput("chip_ctrl_samples",
                       "对照组样本:",
                       choices = sample_names,
                       multiple = TRUE,
                       options = list(create = TRUE))
          ),
          column(6,
            selectizeInput("chip_trt_samples",
                       "处理组样本:",
                       choices = sample_names,
                       multiple = TRUE,
                       options = list(create = TRUE))
          )
        )
      )
    }
  })

  # 显示当前分组状态
  output$chip_current_groups_ui <- renderUI({
    req(chip_data$series_matrix)

    # 检查是否有手动设置的分组
    has_manual <- !is.null(chip_data$manual_ctrl_samples) &&
                  !is.null(chip_data$manual_trt_samples)

    # 检查是否有自动检测的分组
    has_auto <- !is.null(chip_data$group_info) &&
               !is.null(chip_data$group_info$pattern_name)

    if (has_manual) {
      div(
        class = "alert alert-success",
        h5("✅ 当前分组（手动设置）"),
        p(sprintf("对照组 (%d个): %s",
                  length(chip_data$manual_ctrl_samples),
                  paste(chip_data$manual_ctrl_samples, collapse = ", "))),
        p(sprintf("处理组 (%d个): %s",
                  length(chip_data$manual_trt_samples),
                  paste(chip_data$manual_trt_samples, collapse = ", ")))
      )
    } else if (has_auto) {
      div(
        class = "alert alert-info",
        h5("🤖 当前分组（自动检测）"),
        p(sprintf("模式: %s", chip_data$group_info$pattern_name)),
        p(sprintf("对照组 (%d个): %s",
                  length(chip_data$group_info$ctrl_samples),
                  paste(chip_data$group_info$ctrl_samples, collapse = ", "))),
        p(sprintf("处理组 (%d个): %s",
                  length(chip_data$group_info$trt_samples),
                  paste(chip_data$group_info$trt_samples, collapse = ", ")))
      )
    } else {
      div(
        class = "alert alert-warning",
        h5("⚠️ 尚未设置分组"),
        p("请使用上方的快速粘贴功能设置分组，或使用手动选择。")
      )
    }
  })

  # 手动分组UI（折叠状态）
  output$chip_manual_grouping_ui <- renderUI({
    req(chip_data$series_matrix)

    sample_names <- colnames(chip_data$series_matrix)

    tagList(
      h5("📝 手动选择样本（备选方案）", style = "color: #6E6E73;"),
      helpText("如果快速粘贴不方便，可以在这里手动选择样本："),
      fluidRow(
        column(6,
          selectInput("chip_ctrl_samples_manual",
                     "对照组样本:",
                     choices = sample_names,
                     multiple = TRUE)
        ),
        column(6,
          selectInput("chip_trt_samples_manual",
                     "处理组样本:",
                     choices = sample_names,
                     multiple = TRUE)
        )
      ),
      actionButton("chip_apply_manual_groups",
                  "✅ 应用手动选择的分组",
                  class = "btn-success",
                  style = "width: 100%; margin-top: 10px;")
    )
  })

  # 应用手动选择的分组
  observeEvent(input$chip_apply_manual_groups, {
    req(input$chip_ctrl_samples_manual)
    req(input$chip_trt_samples_manual)

    if (length(input$chip_ctrl_samples_manual) == 0 ||
        length(input$chip_trt_samples_manual) == 0) {
      showNotification("请至少为每组选择一个样本！", type = "warning")
      return(NULL)
    }

    chip_data$manual_ctrl_samples <- input$chip_ctrl_samples_manual
    chip_data$manual_trt_samples <- input$chip_trt_samples_manual

    showNotification(
      sprintf("✅ 已应用手动分组: %d 对照 + %d 处理",
              length(input$chip_ctrl_samples_manual),
              length(input$chip_trt_samples_manual)),
      type = "message"
    )
  })

  # 运行差异分析
  observeEvent(input$run_chip_analysis, {
    req(chip_data$series_matrix)

    # 获取分组信息（优先级：手动粘贴 > 自动检测 > 下拉选择）
    if (!is.null(chip_data$manual_ctrl_samples) && !is.null(chip_data$manual_trt_samples)) {
      # 使用手动粘贴的分组
      ctrl_samples <- chip_data$manual_ctrl_samples
      trt_samples <- chip_data$manual_trt_samples
      cat("✅ 使用手动粘贴的分组\n")
    } else if (!is.null(input$chip_use_auto_groups) && input$chip_use_auto_groups &&
               !is.null(chip_data$group_info) && !is.null(chip_data$group_info$pattern_name)) {
      # 使用自动检测的分组
      ctrl_samples <- chip_data$group_info$ctrl_samples
      trt_samples <- chip_data$group_info$trt_samples
      cat("✅ 使用自动检测的分组\n")
    } else {
      # 使用下拉选择的分组
      ctrl_samples <- input$chip_ctrl_samples
      trt_samples <- input$chip_trt_samples
      cat("✅ 使用下拉选择的分组\n")
    }

    # 验证分组
    if (is.null(ctrl_samples) || length(ctrl_samples) == 0 ||
        is.null(trt_samples) || length(trt_samples) == 0) {
      showNotification("请先设置对照组和处理组样本！", type = "error")
      return(NULL)
    }

    # 显示进度
    showNotification("正在运行差异分析...", type = "message")

    # ✅ 优先使用经过完整预处理流程的数据
    if (!is.null(chip_data$standard_expression) && chip_data$ready_for_analysis) {
      # 使用标准格式数据（已探针注释、预处理、去重）
      expr_matrix <- chip_data$standard_expression
      cat("✅ 使用标准格式数据（已探针注释和去重）\n")
      cat(sprintf("   表达矩阵: %d 基因 × %d 样本\n",
                  nrow(expr_matrix), ncol(expr_matrix)))

    } else if (!is.null(chip_data$expr_deduped)) {
      # 使用去重后的数据
      expr_matrix <- chip_data$expr_deduped
      cat("✅ 使用去重后的表达数据\n")

    } else if (!is.null(chip_data$merged_matrix)) {
      # 使用合并后的数据（探针已注释，但未去重）
      # 需要去掉ProbeID和Gene列
      merged_df <- chip_data$merged_matrix
      # 识别数值列（样本数据）
      numeric_cols <- sapply(merged_df, function(x) is.numeric(x))
      expr_matrix <- as.matrix(merged_df[, numeric_cols, drop = FALSE])
      # 使用ProbeID作为行名（因为可能还需要探针级别的信息）
      if ("ProbeID" %in% colnames(merged_df)) {
        rownames(expr_matrix) <- merged_df$ProbeID
      } else {
        rownames(expr_matrix) <- merged_df$Gene
      }
      cat("✅ 使用合并后的数据（探针已注释）\n")

    } else {
      # 最后的备选方案：使用原始Series Matrix
      probe_mapping <- chip_data$probe_mapping

      if (is.null(probe_mapping)) {
        # 如果没有加载注释文件，使用探针ID作为基因符号
        expr_matrix <- chip_data$series_matrix
        cat("⚠️ 未加载注释文件，使用探针ID作为基因符号\n")
      } else {
        # 使用注释映射（旧的自动检测方法）
        expr_matrix <- aggregate_probe_expression(
          chip_data$series_matrix,
          probe_mapping
        )

        if (is.null(expr_matrix)) {
          showNotification("探针注释失败！", type = "error")
          return(NULL)
        }
      }
    }

    # 运行 limma 分析
    limma_res <- run_limma_analysis(
      expr_matrix = expr_matrix,
      ctrl_samples = ctrl_samples,
      trt_samples = trt_samples,
      logfc_threshold = input$chip_logfc_threshold,
      pvalue_threshold = input$chip_pvalue_threshold,
      pval_type = input$chip_pval_type  # 🔧 传入用户选择的P值类型
    )

    if (is.null(limma_res)) {
      showNotification("差异分析失败！", type = "error")
      return(NULL)
    }

    # 转换为标准格式
    formatted_results <- format_chip_results_for_pipeline(
      limma_res,
      expr_matrix,
      ctrl_samples,
      trt_samples
    )

    # 保存结果到 reactiveValues
    chip_data$limma_results <- limma_res
    chip_data$formatted_results <- formatted_results

    showNotification(
      sprintf("✅ 分析完成: %d 个显著差异基因",
              limma_res$n_significant),
      type = "message"
    )
  })

  # 显示结果
  output$chip_results_ui <- renderUI({
    req(chip_data$limma_results)

    limma_res <- chip_data$limma_results

    tagList(
      # 统计摘要
      fluidRow(
        column(3,
          wellPanel(
            style = "background: #f8f9fa; border: 2px solid #dee2e6; text-align: center; padding: 20px;",
            h3(limma_res$n_total, style = "color: #495057; margin: 10px 0;"),
            h6("总基因数", style = "color: #6c757d; margin: 0;")
          )
        ),
        column(3,
          wellPanel(
            style = "background: #e7f3ff; border: 2px solid #007bff; text-align: center; padding: 20px;",
            h3(limma_res$n_significant, style = "color: #007bff; margin: 10px 0;"),
            icon("star", style = "color: #007bff; font-size: 24px;"),
            h6("显著差异", style = "color: #007bff; margin: 5px 0 0 0;")
          )
        ),
        column(3,
          wellPanel(
            style = "background: #d4edda; border: 2px solid #28a745; text-align: center; padding: 20px;",
            h3(limma_res$n_up, style = "color: #28a745; margin: 10px 0;"),
            icon("arrow-up", style = "color: #28a745; font-size: 24px;"),
            h6("上调", style = "color: #28a745; margin: 5px 0 0 0;")
          )
        ),
        column(3,
          wellPanel(
            style = "background: #f8d7da; border: 2px solid #dc3545; text-align: center; padding: 20px;",
            h3(limma_res$n_down, style = "color: #dc3545; margin: 10px 0;"),
            icon("arrow-down", style = "color: #dc3545; font-size: 24px;"),
            h6("下调", style = "color: #dc3545; margin: 5px 0 0 0;")
          )
        )
      ),

      hr(),

      # 结果表格
      h4("差异分析结果"),
      DTOutput("chip_results_table"),

      br(),

      # 下载按钮
      downloadButton("download_chip_results", "📥 下载结果", class = "btn-success")
    )
  })

  # 结果表格
  output$chip_results_table <- renderDT({
    req(chip_data$limma_results)

    results <- chip_data$limma_results$results

    # 🔧 添加显著性标记 - 根据用户选择的P值类型
    pval_col <- if (input$chip_pval_type == "adj.P.Val") "adj.P.Val" else "P.Value"

    results$Significant <- ifelse(
      results[[pval_col]] < input$chip_pvalue_threshold &
        abs(results$logFC) >= input$chip_logfc_threshold,
      "Yes", "No"
    )

    # 🔧 修复：重新排列列顺序，确保ID和SYMBOL列在前
    # 期望顺序：ID, SYMBOL, logFC, AveExpr, t, P.Value, adj.P.Val, B, Significant
    results_ordered <- results[, c("ID", "SYMBOL", "logFC", "AveExpr", "t",
                                    "P.Value", "adj.P.Val", "B", "Significant")]

    # 🔍 诊断：检查ID列的内容
    id_sample <- head(results_ordered$ID, 5)
    cat(sprintf("📊 差异分析结果ID列示例: %s\n", paste(id_sample, collapse = ", ")))
    cat(sprintf("📊 ID列类型: %s\n", class(results_ordered$ID)[1]))

    # 检查ID列是否为基因符号（包含字母）还是EntrezID（纯数字）
    is_entrez_id <- all(grepl("^[0-9]+$", results_ordered$ID[!is.na(results_ordered$ID)]))
    if (!is_entrez_id) {
      cat("⚠️ 警告: ID列包含基因符号而非Entrez Gene ID！\n")
      cat("💡 这可能是因为SOFT文件中缺少EntrezID列\n")
    }

    datatable(
      results_ordered,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(6, 'asc'))  # 按P.Value排序（第6列）
      ),
      filter = 'top',
      rownames = FALSE
    ) %>%
      formatRound(columns = c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val'), digits = 4) %>%
      formatStyle(
        columns = c('ID', 'SYMBOL'),
        backgroundColor = '#e8f4f8',
        fontWeight = 'bold'
      ) %>%
      formatStyle(
        'Significant',
        color = styleEqual(c('Yes', 'No'), c('green', 'grey')),
        fontWeight = 'bold'
      )
  })

  # 下载结果
  output$download_chip_results <- downloadHandler(
    filename = function() {
      sprintf("chip_analysis_results_%s.csv",
              Sys.Date())
    },
    content = function(file) {
      req(chip_data$limma_results)

      results <- chip_data$limma_results$results
      write.csv(results, file, row.names = FALSE)
    }
  )
}

# =====================================================
# 模块结束
# =====================================================

# =====================================================
# 7. 辅助函数：解析粘贴的样本列表
# =====================================================

#' 解析粘贴的样本列表（每行一个样本名）
#'
#' @param pasted_text 粘贴的文本内容
#' @param expr_matrix 表达矩阵（用于验证样本名）
#' @return vector 样本名向量
parse_sample_list <- function(pasted_text, expr_matrix) {
  cat("🔍 开始解析样本列表...\n")

  # 分割成行
  lines <- strsplit(pasted_text, "\n")[[1]]
  lines <- trimws(lines)  # 移除首尾空白
  lines <- lines[lines != ""]  # 移除空行

  if (length(lines) == 0) {
    cat("⚠️  粘贴内容为空\n")
    return(NULL)
  }

  cat(sprintf("📊 读取到 %d 行\n", length(lines)))

  # 获取表达矩阵的样本名
  matrix_samples <- colnames(expr_matrix)

  # 匹配样本名
  valid_samples <- intersect(lines, matrix_samples)

  if (length(valid_samples) == 0) {
    cat("⚠️  未找到匹配的样本\n")
    cat(sprintf("   粘贴的样本: %s\n", paste(head(lines, 3), collapse = ", ")))
    cat(sprintf("   可用样本: %s\n", paste(head(matrix_samples, 3), collapse = ", ")))
    return(NULL)
  }

  cat(sprintf("✅ 匹配成功: %d / %d 个样本\n",
              length(valid_samples), length(lines)))

  # 检查是否有未匹配的样本
  unmatched <- setdiff(lines, matrix_samples)
  if (length(unmatched) > 0) {
    cat(sprintf("⚠️  %d 个样本未匹配: %s\n",
                length(unmatched),
                paste(head(unmatched, 3), collapse = ", ")))
  }

  return(valid_samples)
}

#' 解析用户粘贴的分组信息（旧版本，保留备用）
#'
#' @param pasted_text 粘贴的文本内容
#' @param group_col_name 分组列名（可选）
#' @param expr_matrix 表达矩阵（用于验证样本名）
#' @return list 包含 ctrl_samples 和 trt_samples
parse_pasted_groups <- function(pasted_text, group_col_name, expr_matrix) {
  cat("🔍 开始解析粘贴的分组信息...\n")

  # 分割成行
  lines <- strsplit(pasted_text, "\n")[[1]]
  lines <- lines[lines != ""]  # 移除空行

  if (length(lines) == 0) {
    return(list(success = FALSE, error = "粘贴内容为空"))
  }

  # 方法1: 如果指定了分组列名，按列名解析
  if (!is.null(group_col_name) && group_col_name != "") {
    cat(sprintf("📋 使用分组列名: %s\n", group_col_name))

    # 尝试读取为表格
    tryCatch({
      # 读取文本
      text_conn <- textConnection(lines)
      df <- read.table(text_conn, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE,
                       quote = "\"", comment.char = "")
      close(text_conn)

      # 检查是否存在分组列
      if (!group_col_name %in% colnames(df)) {
        return(list(success = FALSE,
                    error = sprintf("未找到分组列 '%s'，可用列: %s",
                                   group_col_name,
                                   paste(colnames(df), collapse = ", "))))
      }

      # 提取分组信息
      groups <- df[[group_col_name]]

      # 获取唯一分组
      unique_groups <- unique(groups)

      if (length(unique_groups) < 2) {
        return(list(success = FALSE, error = "分组列中只有一个组，需要至少2个组"))
      }

      # 如果超过2个组，取前2个
      if (length(unique_groups) > 2) {
        cat(sprintf("⚠️  检测到 %d 个组，将使用前两个: %s\n",
                    length(unique_groups),
                    paste(unique_groups[1:2], collapse = ", ")))
        unique_groups <- unique_groups[1:2]
      }

      # 根据分组列提取样本
      ctrl_name <- unique_groups[1]
      trt_name <- unique_groups[2]

      ctrl_samples <- colnames(df)[groups == ctrl_name]
      trt_samples <- colnames(df)[groups == trt_name]

      # 移除分组列本身
      ctrl_samples <- ctrl_samples[ctrl_samples != group_col_name]
      trt_samples <- trt_samples[trt_samples != group_col_name]

      cat(sprintf("✅ 解析成功: %d 对照 (%s) + %d 处理 (%s)\n",
                  length(ctrl_samples), ctrl_name,
                  length(trt_samples), trt_name))

      return(list(
        success = TRUE,
        ctrl_samples = ctrl_samples,
        trt_samples = trt_samples
      ))

    }, error = function(e) {
      return(list(success = FALSE, error = paste("解析表格失败:", e$message)))
    })
  }

  # 方法2: 没有指定列名，自动检测
  cat("🤖 自动检测分组模式...\n")

  # 尝试检测第一行是否为列名
  first_line <- lines[1]

  tryCatch({
    # 读取为表格
    text_conn <- textConnection(lines)
    df <- read.table(text_conn, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "\"", comment.char = "")
    close(text_conn)

    # 获取所有样本名（列名，排除第一列ID列）
    all_samples <- colnames(df)[-1]  # 排除ID列

    # 获取表达矩阵的样本名
    matrix_samples <- colnames(expr_matrix)

    # 找到交集
    common_samples <- intersect(all_samples, matrix_samples)

    if (length(common_samples) < 2) {
      return(list(success = FALSE,
                  error = sprintf("在粘贴内容中只找到 %d 个有效样本，需要至少2个",
                                 length(common_samples))))
    }

    cat(sprintf("📊 找到 %d 个有效样本\n", length(common_samples)))

    # 尝试从列名中自动分组
    # 使用现有的自动检测函数
    sample_names <- common_samples

    # 如果有多余的列，尝试从第二列开始检测分组模式
    if (ncol(df) > 2) {
      # 检查第二列是否有分组信息
      potential_groups <- df[[2]]

      if (length(unique(potential_groups)) == 2) {
        unique_g <- unique(potential_groups)
        ctrl_samples <- sample_names[potential_groups == unique_g[1]]
        trt_samples <- sample_names[potential_groups == unique_g[2]]

        cat(sprintf("✅ 自动检测分组: %d 对照 + %d 处理\n",
                    length(ctrl_samples), length(trt_samples)))

        return(list(
          success = TRUE,
          ctrl_samples = ctrl_samples,
          trt_samples = trt_samples
        ))
      }
    }

    # 如果无法自动检测，返回所有样本让用户手动选择
    cat("⚠️  无法自动检测分组，返回所有样本供手动选择\n")

    # 默认：前一半作为对照，后一半作为处理
    n_samples <- length(sample_names)
    mid_point <- ceiling(n_samples / 2)

    ctrl_samples <- sample_names[1:mid_point]
    trt_samples <- sample_names[(mid_point+1):n_samples]

    cat(sprintf("⚠️  默认分组: 前 %d 个为对照，后 %d 个为处理\n",
                length(ctrl_samples), length(trt_samples)))

    return(list(
      success = TRUE,
      ctrl_samples = ctrl_samples,
      trt_samples = trt_samples,
      warning = "无法自动检测分组模式，已按位置默认分组，请手动调整"
    ))

  }, error = function(e) {
    return(list(success = FALSE, error = paste("解析失败:", e$message)))
  })
}

