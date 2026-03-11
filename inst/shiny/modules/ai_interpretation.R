# =====================================================
# AI 解读模块
# 支持多种 API 提供商（OpenAI、智谱AI、本地模型等）
# =====================================================

library(httr)
library(jsonlite)
library(base64enc)

# =====================================================
# API 提供商配置
# =====================================================

api_providers <- list(
  openai = list(
    name = "OpenAI",
    endpoint = "https://api.openai.com/v1/chat/completions",
    models = c("gpt-4o", "gpt-4-turbo", "gpt-4", "gpt-3.5-turbo"),
    default_model = "gpt-4o"
  ),
  zhipu = list(
    name = "智谱AI (GLM)",
    endpoint = "https://open.bigmodel.cn/api/paas/v4/chat/completions",
    models = c("glm-4", "glm-4-flash", "glm-4-plus"),
    default_model = "glm-4-flash"
  ),
  deepseek = list(
    name = "DeepSeek",
    endpoint = "https://api.deepseek.com/v1/chat/completions",
    models = c("deepseek-chat", "deepseek-coder"),
    default_model = "deepseek-chat"
  ),
  local = list(
    name = "本地模型",
    endpoint = "http://localhost:8000/v1/chat/completions",
    models = c("custom"),
    default_model = "custom"
  ),
  custom = list(
    name = "自定义API",
    endpoint = "",
    models = c("custom"),
    default_model = "custom"
  )
)

# =====================================================
# 数据汇总函数
# =====================================================

#' 汇总差异分析结果
summarize_deg_results <- function(deg_data) {
  if (is.null(deg_data) || is.null(deg_data$deg_df)) {
    return(NULL)
  }

  deg_df <- deg_data$deg_df

  # 基本统计
  stats <- list(
    total_genes = nrow(deg_df),
    up_genes = sum(deg_df$Status == "Up", na.rm = TRUE),
    down_genes = sum(deg_df$Status == "Down", na.rm = TRUE),
    significant_genes = sum(deg_df$Status != "Not DE", na.rm = TRUE)
  )

  # Top 上调基因
  top_up <- deg_df %>%
    dplyr::filter(Status == "Up") %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    head(15) %>%
    dplyr::select(SYMBOL, log2FoldChange, padj)

  # Top 下调基因
  top_down <- deg_df %>%
    dplyr::filter(Status == "Down") %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    head(15) %>%
    dplyr::select(SYMBOL, log2FoldChange, padj)

  # log2FC 分布统计
  log2fc_stats <- list(
    min = round(min(deg_df$log2FoldChange, na.rm = TRUE), 3),
    max = round(max(deg_df$log2FoldChange, na.rm = TRUE), 3),
    median = round(median(deg_df$log2FoldChange, na.rm = TRUE), 3),
    mean = round(mean(deg_df$log2FoldChange, na.rm = TRUE), 3)
  )

  return(list(
    stats = stats,
    top_up = top_up,
    top_down = top_down,
    log2fc_stats = log2fc_stats
  ))
}

#' 汇总富集分析结果
summarize_enrichment_results <- function(enrichment_df, analysis_type = "KEGG") {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }

  # 确保数据框格式正确
  if (!is.data.frame(enrichment_df)) {
    # 尝试从 reactive 中提取
    tryCatch({
      enrichment_df <- enrichment_df()
    }, error = function(e) {
      return(NULL)
    })
  }

  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }

  # Top 富集通路
  top_pathways <- enrichment_df %>%
    dplyr::arrange(p.adjust) %>%
    head(15)

  # 选择可用列
  available_cols <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                                 "pvalue", "p.adjust", "geneID", "Count"),
                               colnames(top_pathways))
  top_pathways <- top_pathways[, available_cols, drop = FALSE]

  # 统计信息
  stats <- list(
    total_pathways = nrow(enrichment_df),
    significant_pathways = sum(enrichment_df$p.adjust < 0.05, na.rm = TRUE),
    top_pathway = if(nrow(top_pathways) > 0) top_pathways$Description[1] else "无",
    min_pvalue = if(nrow(enrichment_df) > 0) round(min(enrichment_df$pvalue, na.rm = TRUE), 4) else NA
  )

  return(list(
    type = analysis_type,
    stats = stats,
    top_pathways = top_pathways
  ))
}

#' 汇总 TF/通路活性结果
summarize_activity_results <- function(activity_df, analysis_type = "TF") {
  if (is.null(activity_df)) {
    return(NULL)
  }

  # 尝试从 reactive 中提取
  if (!is.data.frame(activity_df)) {
    tryCatch({
      activity_df <- activity_df()
    }, error = function(e) {
      return(NULL)
    })
  }

  if (is.null(activity_df) || nrow(activity_df) == 0) {
    return(NULL)
  }

  # 确定分数列名
  score_col <- intersect(c("score", "activity", "estimate"), colnames(activity_df))[1]
  pval_col <- intersect(c("p_value", "p.value", "pval", "p.adjust"), colnames(activity_df))[1]
  source_col <- intersect(c("source", "tf", "pathway", "ID"), colnames(activity_df))[1]

  if (is.na(score_col) || is.na(source_col)) {
    return(NULL)
  }

  # Top 激活/抑制
  top_activated <- activity_df %>%
    dplyr::filter(.data[[score_col]] > 0) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]])) %>%
    head(10) %>%
    dplyr::select(dplyr::all_of(c(source_col, score_col, pval_col)))

  top_inhibited <- activity_df %>%
    dplyr::filter(.data[[score_col]] < 0) %>%
    dplyr::arrange(.data[[score_col]]) %>%
    head(10) %>%
    dplyr::select(dplyr::all_of(c(source_col, score_col, pval_col)))

  stats <- list(
    total = nrow(activity_df),
    activated = sum(activity_df[[score_col]] > 0, na.rm = TRUE),
    inhibited = sum(activity_df[[score_col]] < 0, na.rm = TRUE),
    significant = if (!is.na(pval_col)) sum(activity_df[[pval_col]] < 0.05, na.rm = TRUE) else NA
  )

  return(list(
    type = analysis_type,
    stats = stats,
    top_activated = top_activated,
    top_inhibited = top_inhibited
  ))
}

# =====================================================
# Prompt 生成函数
# =====================================================

#' 生成 AI 解读 Prompt
generate_ai_prompt <- function(deg_summary, kegg_summary, go_summary,
                                tf_summary, pathway_summary,
                                analysis_type = "comprehensive",
                                language = "中文",
                                sample_info = NULL) {

  lang_instruction <- if (language == "中文") {
    "请使用中文回答，使用专业的生物信息学术语。"
  } else {
    "Please respond in English using professional bioinformatics terminology."
  }

  # 构建基础 prompt
  prompt_parts <- list()

  # 样本背景信息
  sample_context <- ""
  if (!is.null(sample_info)) {
    sample_parts <- c()
    if (!is.null(sample_info$organism) && sample_info$organism != "") {
      sample_parts <- c(sample_parts, paste0("- 物种: ", sample_info$organism))
    }
    if (!is.null(sample_info$tissue) && sample_info$tissue != "") {
      sample_parts <- c(sample_parts, paste0("- 组织/细胞: ", sample_info$tissue))
    }
    if (!is.null(sample_info$condition) && sample_info$condition != "") {
      sample_parts <- c(sample_parts, paste0("- 实验处理: ", sample_info$condition))
    }
    if (!is.null(sample_info$control) && sample_info$control != "") {
      sample_parts <- c(sample_parts, paste0("- 对照组: ", sample_info$control))
    }
    if (!is.null(sample_info$description) && sample_info$description != "") {
      sample_parts <- c(sample_parts, paste0("- 补充说明: ", sample_info$description))
    }

    if (length(sample_parts) > 0) {
      sample_context <- paste0(
        "## 样本背景信息\n\n",
        paste(sample_parts, collapse = "\n"),
        "\n\n**请结合上述样本背景信息进行针对性解读，重点关注与实验处理相关的生物学变化。**\n\n"
      )
    }
  }

  prompt_parts[[1]] <- paste0(
    "你是一位专业的生物信息学分析师，请基于以下RNA-seq差异表达分析结果，",
    "生成一份专业、全面的生物学解读报告。\n\n",
    lang_instruction, "\n\n",
    "---\n\n",
    "# 分析数据\n\n",
    sample_context
  )

  # 差异分析数据
  if (!is.null(deg_summary) && analysis_type %in% c("comprehensive", "deg_only")) {
    stats <- deg_summary$stats
    prompt_parts[[length(prompt_parts) + 1]] <- paste0(
      "## 1. 差异表达分析结果\n\n",
      "- 总检测基因数: ", stats$total_genes, "\n",
      "- 显著差异基因数: ", stats$significant_genes, "\n",
      "  - 上调基因: ", stats$up_genes, "\n",
      "  - 下调基因: ", stats$down_genes, "\n",
      "- log2FC 范围: [", deg_summary$log2fc_stats$min, ", ", deg_summary$log2fc_stats$max, "]\n",
      "- log2FC 中位数: ", deg_summary$log2fc_stats$median, "\n\n",
      "### Top 15 上调基因\n",
      "```\n",
      paste(capture.output(print(deg_summary$top_up, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n",
      "### Top 15 下调基因\n",
      "```\n",
      paste(capture.output(print(deg_summary$top_down, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n"
    )
  }

  # KEGG 富集数据
  if (!is.null(kegg_summary) && analysis_type %in% c("comprehensive", "enrichment_only")) {
    stats <- kegg_summary$stats
    prompt_parts[[length(prompt_parts) + 1]] <- paste0(
      "## 2. KEGG 富集分析结果\n\n",
      "- 富集通路总数: ", stats$total_pathways, "\n",
      "- 显著通路数 (p.adj < 0.05): ", stats$significant_pathways, "\n\n",
      "### Top 15 富集通路\n",
      "```\n",
      paste(capture.output(print(kegg_summary$top_pathways, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n"
    )
  }

  # GO 富集数据
  if (!is.null(go_summary) && analysis_type %in% c("comprehensive", "enrichment_only")) {
    stats <- go_summary$stats
    prompt_parts[[length(prompt_parts) + 1]] <- paste0(
      "## 3. GO 富集分析结果\n\n",
      "- GO Terms 总数: ", stats$total_pathways, "\n",
      "- 显著 GO Terms 数: ", stats$significant_pathways, "\n\n",
      "### Top 15 GO Terms\n",
      "```\n",
      paste(capture.output(print(go_summary$top_pathways, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n"
    )
  }

  # TF 活性数据
  if (!is.null(tf_summary) && analysis_type %in% c("comprehensive", "tf_only")) {
    stats <- tf_summary$stats
    prompt_parts[[length(prompt_parts) + 1]] <- paste0(
      "## 4. 转录因子活性分析\n\n",
      "- 分析的 TF 总数: ", stats$total, "\n",
      "- 激活的 TF 数: ", stats$activated, "\n",
      "- 抑制的 TF 数: ", stats$inhibited, "\n\n",
      "### Top 10 激活的转录因子\n",
      "```\n",
      paste(capture.output(print(tf_summary$top_activated, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n",
      "### Top 10 抑制的转录因子\n",
      "```\n",
      paste(capture.output(print(tf_summary$top_inhibited, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n"
    )
  }

  # 通路活性数据
  if (!is.null(pathway_summary) && analysis_type %in% c("comprehensive", "pathway_only")) {
    stats <- pathway_summary$stats
    prompt_parts[[length(prompt_parts) + 1]] <- paste0(
      "## 5. 通路活性分析\n\n",
      "- 分析的通路总数: ", stats$total, "\n",
      "- 激活的通路数: ", stats$activated, "\n",
      "- 抑制的通路数: ", stats$inhibited, "\n\n",
      "### Top 10 激活的通路\n",
      "```\n",
      paste(capture.output(print(pathway_summary$top_activated, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n",
      "### Top 10 抑制的通路\n",
      "```\n",
      paste(capture.output(print(pathway_summary$top_inhibited, row.names = FALSE)), collapse = "\n"),
      "\n```\n\n"
    )
  }

  # 输出要求
  prompt_parts[[length(prompt_parts) + 1]] <- paste0(
    "---\n\n",
    "# 输出要求\n\n",
    "请按照以下结构生成解读报告：\n\n",
    "## 1. 关键发现总结\n",
    "简明扼要地总结3-5个最重要的生物学发现。\n\n",
    "## 2. 差异表达模式分析\n",
    "分析上调和下调基因的整体趋势，讨论Top基因的潜在功能意义。\n\n",
    "## 3. 通路富集解读\n",
    "解读显著富集的通路，分析通路之间的关联性，讨论与实验处理的潜在关系。\n\n",
    "## 4. 转录调控机制（如有数据）\n",
    "分析活性变化的转录因子，讨论TF与靶基因的一致性。\n\n",
    "## 5. 生物学意义\n",
    "整合所有分析结果，提出可能的分子机制假设。\n\n",
    "## 6. 后续研究建议\n",
    "建议需要验证的关键基因/通路，推荐实验验证方法。\n\n",
    "---\n\n",
    "**注意**: 请使用Markdown格式输出，适当引用具体数据和p值，客观分析，避免过度解读。"
  )

  return(paste(prompt_parts, collapse = ""))
}

# =====================================================
# API 调用函数
# =====================================================

#' 调用 AI API
call_ai_api <- function(prompt, config) {

  # 验证配置
  if (is.null(config$api_key) || config$api_key == "") {
    return(list(
      success = FALSE,
      error = "API Key 未配置，请先设置 API Key"
    ))
  }

  provider <- config$provider
  provider_config <- api_providers[[provider]]

  if (is.null(provider_config)) {
    return(list(
      success = FALSE,
      error = paste("未知的 API 提供商:", provider)
    ))
  }

  # 确定端点
  endpoint <- config$api_endpoint
  if (is.null(endpoint) || endpoint == "") {
    endpoint <- provider_config$endpoint
  }

  tryCatch({
    # 构建请求体（OpenAI 兼容格式）
    request_body <- list(
      model = config$model,
      messages = list(
        list(role = "system", content = "你是一位专业的生物信息学分析师，擅长RNA-seq数据分析和生物学解读。"),
        list(role = "user", content = prompt)
      ),
      max_tokens = config$max_tokens %||% 4000,
      temperature = config$temperature %||% 0.7
    )

    # 发送请求
    response <- POST(
      url = endpoint,
      add_headers(
        `Content-Type` = "application/json",
        `Authorization` = paste("Bearer", config$api_key)
      ),
      body = toJSON(request_body, auto_unbox = TRUE),
      encode = "raw",
      timeout(120)
    )

    # 检查响应状态
    if (http_error(response)) {
      error_content <- content(response, "text", encoding = "UTF-8")
      return(list(
        success = FALSE,
        error = paste("API 请求失败:", status_code(response), error_content)
      ))
    }

    # 解析响应 - 使用 simplifyVector = FALSE 保持列表结构，避免 data.frame 访问问题
    result <- fromJSON(content(response, "text", encoding = "UTF-8"), simplifyVector = FALSE)

    # 提取内容（兼容不同 API 格式）
    # 注意: fromJSON 默认会将 JSON 数组转换为 data.frame
    # 使用 simplifyVector = FALSE 后，result$choices 是一个列表

    content_text <- NULL
    tokens_used <- NA

    # 方式1: OpenAI/DeepSeek/智谱AI 标准格式 (choices 是列表)
    if (!is.null(result$choices) && length(result$choices) > 0) {
      choice <- result$choices[[1]]  # 现在是列表，[[1]] 正确获取第一个元素

      # 检查 message 结构
      if (!is.null(choice$message)) {
        content_text <- choice$message$content
      } else if (!is.null(choice$text)) {
        # 某些旧版 API 可能直接返回 text
        content_text <- choice$text
      }

      # 获取 token 使用量
      if (!is.null(result$usage)) {
        tokens_used <- if (!is.null(result$usage$total_tokens)) {
          result$usage$total_tokens
        } else {
          NA
        }
      }
    }
    # 方式2: 智谱AI 旧格式 (data 包装)
    else if (!is.null(result$data)) {
      if (!is.null(result$data$choices) && length(result$data$choices) > 0) {
        choice <- result$data$choices[[1]]
        if (!is.null(choice$message)) {
          content_text <- choice$message$content
        } else if (!is.null(choice$content)) {
          content_text <- choice$content
        }
      }
      if (!is.null(result$data$usage)) {
        tokens_used <- result$data$usage$total_tokens
      }
    }
    # 方式3: 错误响应
    else if (!is.null(result$error)) {
      error_msg <- if (is.list(result$error)) {
        result$error$message %||% toString(result$error)
      } else {
        toString(result$error)
      }
      return(list(
        success = FALSE,
        error = paste("API 返回错误:", error_msg)
      ))
    }

    # 检查是否成功获取内容
    if (is.null(content_text) || content_text == "") {
      return(list(
        success = FALSE,
        error = paste("无法解析 API 响应结构。响应预览:",
                      substr(toString(result), 1, 500))
      ))
    }

    return(list(
      success = TRUE,
      content = content_text,
      tokens_used = tokens_used,
      model = config$model
    ))

  }, error = function(e) {
    return(list(
      success = FALSE,
      error = paste("API 调用出错:", e$message)
    ))
  })
}

# =====================================================
# API 配置文件路径（用户配置）
# =====================================================

get_api_config_file <- function() {
  # 配置文件位置: 用户主目录/.yuanseq/api_config.json
  # Windows: C:\Users\{用户名}\.yuanseq\api_config.json
  # Mac/Linux: ~/.yuanseq/api_config.json
  file.path(path.expand("~"), ".yuanseq", "api_config.json")
}

load_api_config_from_file <- function() {
  config_file <- get_api_config_file()
  if (!file.exists(config_file)) {
    return(NULL)
  }
  tryCatch({
    config <- jsonlite::fromJSON(config_file)
    return(config)
  }, error = function(e) {
    return(NULL)
  })
}

# =====================================================
# AI 解读 Server 模块
# =====================================================

ai_interpretation_server <- function(input, output, session,
                                      deg_results,
                                      kegg_results = NULL,
                                      go_results = NULL,
                                      tf_activity_results = NULL,
                                      pathway_activity_results = NULL,
                                      # 图片文件路径（用于导出HTML时嵌入）
                                      volcano_plot_path = NULL,
                                      heatmap_plot_path = NULL,
                                      kegg_plot_path = NULL,
                                      go_plot_path = NULL,
                                      tf_plot_path = NULL,
                                      pathway_plot_path = NULL) {

  # 历史记录存储
  interpretation_history <- reactiveVal(data.frame(
    timestamp = character(),
    type = character(),
    model = character(),
    tokens = integer(),
    stringsAsFactors = FALSE
  ))

  # 存储图片路径
  plot_paths <- reactiveValues(
    volcano = NULL,
    heatmap = NULL,
    kegg = NULL,
    go = NULL,
    tf = NULL,
    pathway = NULL
  )

  # ========== 从配置文件加载 API 设置 ==========

  saved_config <- reactiveVal(NULL)

  # 在session初始化时加载配置
  observe({
    config <- load_api_config_from_file()
    if (!is.null(config)) {
      saved_config(config)
      # 更新输入框
      if (!is.null(config$provider)) {
        updateSelectInput(session, "ai_provider", selected = config$provider)
      }
      if (!is.null(config$api_key)) {
        updateTextInput(session, "ai_api_key", value = config$api_key)
      }
      if (!is.null(config$custom_endpoint)) {
        updateTextInput(session, "ai_custom_endpoint", value = config$custom_endpoint)
      }
      if (!is.null(config$model)) {
        # 延迟更新模型选择（等待UI渲染）
        delay(100, {
          updateSelectInput(session, "ai_model", selected = config$model)
        })
      }
    }
  }, priority = 1000)

  # 配置状态显示
  output$api_config_status <- renderUI({
    config <- saved_config()
    config_file <- get_api_config_file()

    if (!is.null(config) && !is.null(config$api_key) && config$api_key != "" && !grepl("在这里", config$api_key)) {
      tags$div(
        class = "alert alert-success small",
        style = "margin-top: 10px; margin-bottom: 0; padding: 8px;",
        tags$strong("✅ 已加载配置文件"),
        tags$br(),
        tags$span(class = "text-muted", paste("提供商:", config$provider, "| 模型:", config$model))
      )
    } else {
      tags$div(
        class = "alert alert-warning small",
        style = "margin-top: 10px; margin-bottom: 0; padding: 8px;",
        tags$strong("⚠️ 未检测到有效配置"),
        tags$br(),
        tags$span("请将 ", tags$code("api_config.example.json"), " 复制为 "),
        tags$code("api_config.json"),
        tags$br(),
        tags$span("填入API Key后保存到: "),
        tags$code(config_file)
      )
    }
  })

  # API 配置 reactive
  api_config <- reactive({
    # 优先使用配置文件，其次使用界面输入
    saved <- saved_config()

    # 获取 API Key
    api_key <- input$ai_api_key
    if (is.null(api_key) || api_key == "") {
      if (!is.null(saved$api_key) && !grepl("在这里", saved$api_key)) {
        api_key <- saved$api_key
      } else {
        api_key <- Sys.getenv("YUANSEQ_AI_API_KEY", unset = "")
      }
    }

    # 获取提供商
    provider <- input$ai_provider
    if (is.null(provider) || provider == "") {
      provider <- saved$provider %||% "deepseek"
    }

    # 获取模型
    model <- input$ai_model
    if (is.null(model) || model == "") {
      model <- saved$model %||% "deepseek-chat"
    }

    # 获取端点
    endpoint <- NULL
    if (provider == "custom") {
      endpoint <- input$ai_custom_endpoint
      if (is.null(endpoint) || endpoint == "") {
        endpoint <- saved$custom_endpoint
      }
    }

    list(
      provider = provider,
      api_key = api_key,
      api_endpoint = endpoint,
      model = model,
      max_tokens = input$ai_max_tokens %||% 4000,
      temperature = input$ai_temperature %||% 0.7
    )
  })

  # 动态模型选择器
  output$ai_model_selector <- renderUI({
    provider <- input$ai_provider
    models <- api_providers[[provider]]$models
    default_model <- api_providers[[provider]]$default_model

    selectInput("ai_model", "模型选择",
                choices = models,
                selected = default_model)
  })

  # 数据状态显示
  output$ai_data_summary <- renderUI({
    # 检查各数据源是否可用
    has_deg <- !is.null(deg_results) && !is.null(tryCatch(deg_results(), error = function(e) NULL))
    has_kegg <- !is.null(kegg_results) && !is.null(tryCatch(kegg_results(), error = function(e) NULL))
    has_go <- !is.null(go_results) && !is.null(tryCatch(go_results(), error = function(e) NULL))

    # TF 活性结果需要特殊处理
    has_tf <- FALSE
    if (!is.null(tf_activity_results)) {
      tryCatch({
        tf_data <- tf_activity_results()
        has_tf <- !is.null(tf_data) && nrow(tf_data) > 0
      }, error = function(e) {})
    }

    # 通路活性结果
    has_pathway <- FALSE
    if (!is.null(pathway_activity_results)) {
      tryCatch({
        pathway_data <- pathway_activity_results()
        has_pathway <- !is.null(pathway_data) && nrow(pathway_data) > 0
      }, error = function(e) {})
    }

    tagList(
      h4("📋 数据状态"),

      tags$div(
        class = "alert alert-info",
        style = "padding: 15px;",

        tags$table(
          style = "width: 100%;",

          tags$tr(
            tags$td(style = "padding: 8px; width: 60%;", "📊 差异分析结果"),
            tags$td(style = "padding: 8px;",
                    if(has_deg) {
                      tryCatch({
                        deg_data <- deg_results()
                        tags$span(class = "text-success",
                                  paste0("✅ ", nrow(deg_data$deg_df), " 个基因"))
                      }, error = function(e) {
                        tags$span(class = "text-muted", "❌ 未运行")
                      })
                    } else {
                      tags$span(class = "text-muted", "❌ 未运行")
                    })
          ),

          tags$tr(
            tags$td(style = "padding: 8px;", "🧬 KEGG 富集结果"),
            tags$td(style = "padding: 8px;",
                    if(has_kegg) {
                      tryCatch({
                        kegg_data <- kegg_results()
                        tags$span(class = "text-success",
                                  paste0("✅ ", nrow(kegg_data), " 个通路"))
                      }, error = function(e) {
                        tags$span(class = "text-muted", "❌ 未运行")
                      })
                    } else {
                      tags$span(class = "text-muted", "❌ 未运行")
                    })
          ),

          tags$tr(
            tags$td(style = "padding: 8px;", "🧬 GO 富集结果"),
            tags$td(style = "padding: 8px;",
                    if(has_go) {
                      tryCatch({
                        go_data <- go_results()
                        tags$span(class = "text-success",
                                  paste0("✅ ", nrow(go_data), " 个 GO term"))
                      }, error = function(e) {
                        tags$span(class = "text-muted", "❌ 未运行")
                      })
                    } else {
                      tags$span(class = "text-muted", "❌ 未运行")
                    })
          ),

          tags$tr(
            tags$td(style = "padding: 8px;", "🔬 TF 活性分析"),
            tags$td(style = "padding: 8px;",
                    if(has_tf) {
                      tags$span(class = "text-success", "✅ 已完成")
                    } else {
                      tags$span(class = "text-muted", "❌ 未运行")
                    })
          ),

          tags$tr(
            tags$td(style = "padding: 8px;", "🛤️ 通路活性分析"),
            tags$td(style = "padding: 8px;",
                    if(has_pathway) {
                      tags$span(class = "text-success", "✅ 已完成")
                    } else {
                      tags$span(class = "text-muted", "❌ 未运行")
                    })
          )
        )
      ),

      if (!has_deg) {
        tags$div(
          class = "alert alert-warning",
          tags$strong("⚠️ 提示: "),
          "差异分析是 AI 解读的基础，请先在「数据概览」标签页运行差异分析。"
        )
      } else {
        NULL
      }
    )
  })

  # AI 解读结果
  ai_result <- reactiveVal(NULL)

  # 进度状态
  ai_progress <- reactiveVal(list(
    active = FALSE,
    step = 0,
    total = 5,
    message = ""
  ))

  # 进度轮询 - 每500ms检查一次进度状态
  observe({
    invalidateLater(500, session)
    progress <- ai_progress()
    # 当进度激活时，强制触发 ai_interpretation_result 重新渲染
    if (progress$active) {
      ai_progress(progress)  # 触发依赖更新
    }
  })

  # 更新进度的辅助函数
  update_progress <- function(step, message) {
    ai_progress(list(
      active = TRUE,
      step = step,
      total = 5,
      message = message
    ))
  }

  # 运行 AI 解读
  observeEvent(input$run_ai_interpretation, {

    # 验证数据
    deg_data <- tryCatch(deg_results(), error = function(e) NULL)
    if (is.null(deg_data)) {
      showNotification("❌ 请先运行差异分析", type = "error", duration = 5)
      return()
    }

    # 验证 API Key
    config <- api_config()
    if (is.null(config$api_key) || config$api_key == "") {
      showNotification("❌ 请先配置 API Key", type = "error", duration = 5)
      return()
    }

    # 重置结果并开始进度
    ai_result(NULL)
    update_progress(1, "📊 正在汇总差异分析数据...")

    # 汇总数据
    deg_summary <- summarize_deg_results(deg_data)

    update_progress(2, "🧬 正在汇总富集分析结果...")
    Sys.sleep(0.3)  # 让 UI 有时间更新

    kegg_summary <- NULL
    if (!is.null(kegg_results)) {
      tryCatch({
        kegg_data <- kegg_results()
        kegg_summary <- summarize_enrichment_results(kegg_data, "KEGG")
      }, error = function(e) {})
    }

    go_summary <- NULL
    if (!is.null(go_results)) {
      tryCatch({
        go_data <- go_results()
        go_summary <- summarize_enrichment_results(go_data, "GO")
      }, error = function(e) {})
    }

    update_progress(3, "🔬 正在汇总活性分析结果...")
    Sys.sleep(0.3)

    tf_summary <- NULL
    if (!is.null(tf_activity_results)) {
      tryCatch({
        tf_data <- tf_activity_results()
        tf_summary <- summarize_activity_results(tf_data, "TF")
      }, error = function(e) {})
    }

    pathway_summary <- NULL
    if (!is.null(pathway_activity_results)) {
      tryCatch({
        pathway_data <- pathway_activity_results()
        pathway_summary <- summarize_activity_results(pathway_data, "Pathway")
      }, error = function(e) {})
    }

    update_progress(4, "📝 正在生成分析提示词...")
    Sys.sleep(0.3)

    # 收集样本信息
    sample_info <- list(
      organism = input$sample_organism,
      tissue = input$sample_tissue,
      condition = input$sample_condition,
      control = input$sample_control,
      description = input$sample_description
    )

    # 生成 Prompt
    prompt <- generate_ai_prompt(
      deg_summary = deg_summary,
      kegg_summary = kegg_summary,
      go_summary = go_summary,
      tf_summary = tf_summary,
      pathway_summary = pathway_summary,
      analysis_type = input$ai_analysis_type,
      language = input$ai_language,
      sample_info = sample_info
    )

    update_progress(5, paste0("🤖 正在调用 ", config$model, " API，请耐心等待..."))

    # 调用 API
    result <- call_ai_api(prompt, config)

    # 清除进度
    ai_progress(list(active = FALSE, step = 0, total = 5, message = ""))

    if (result$success) {
      ai_result(result)

      # 更新历史记录
      history <- interpretation_history()
      new_row <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        type = input$ai_analysis_type,
        model = config$model,
        tokens = result$tokens_used %||% NA,
        stringsAsFactors = FALSE
      )
      interpretation_history(rbind(history, new_row))

      showNotification(paste0("✅ AI 解读完成! 使用 Token: ", result$tokens_used %||% "N/A"),
                       type = "message", duration = 5)
    } else {
      showNotification(paste0("❌ AI 解读失败: ", result$error), type = "error", duration = 10)

      # 提供常见错误解决方案
      if (grepl("401|Unauthorized|invalid", result$error, ignore.case = TRUE)) {
        showNotification("💡 提示: API Key 无效或已过期，请检查配置", type = "warning", duration = 10)
      } else if (grepl("429|rate limit", result$error, ignore.case = TRUE)) {
        showNotification("💡 提示: API 调用频率超限，请稍后重试", type = "warning", duration = 10)
      } else if (grepl("timeout|connection", result$error, ignore.case = TRUE)) {
        showNotification("💡 提示: 网络连接失败，请检查网络或使用代理", type = "warning", duration = 10)
      }
    }
  })

  # 显示 AI 解读结果
  output$ai_interpretation_result <- renderUI({

    # 检查进度状态
    progress <- ai_progress()

    if (progress$active) {
      # 显示进度条
      pct <- round(progress$step / progress$total * 100)
      return(
        tags$div(
          class = "alert alert-info",
          style = "padding: 20px;",
          tags$h4("🤖 AI 解读进行中..."),
          tags$p(progress$message),
          tags$div(
            class = "progress",
            style = "height: 25px; margin-top: 15px;",
            tags$div(
              class = "progress-bar progress-bar-striped active",
              role = "progressbar",
              style = paste0("width: ", pct, "%;"),
              `aria-valuenow` = pct,
              `aria-valuemin` = "0",
              `aria-valuemax` = "100",
              paste0(pct, "%")
            )
          ),
          tags$p(class = "text-muted small", style = "margin-top: 10px;",
                 "步骤 ", progress$step, " / ", progress$total,
                 " - API 调用可能需要 10-60 秒，请耐心等待")
        )
      )
    }

    result <- ai_result()

    if (is.null(result)) {
      return(
        tags$div(
          class = "alert alert-secondary",
          style = "text-align: center; padding: 40px;",
          tags$h4("🤖 AI 解读"),
          tags$p("配置好 API Key 后，点击「生成 AI 解读」按钮开始分析"),
          tags$hr(),
          tags$p(class = "text-muted small",
                 "支持 OpenAI、智谱AI、DeepSeek、本地模型等多种 API 服务")
        )
      )
    }

    if (!result$success) {
      return(
        tags$div(
          class = "alert alert-danger",
          tags$h4("❌ 解读失败"),
          tags$p(result$error)
        )
      )
    }

    # 成功结果 - 渲染 Markdown
    tagList(
      tags$div(
        class = "alert alert-success",
        tags$h4("✅ AI 解读结果"),
        tags$p(paste0("模型: ", result$model, " | Token: ", result$tokens_used %||% "N/A"))
      ),

      tags$div(
        class = "ai-result-content",
        style = "background: #fafafa; padding: 20px; border-radius: 8px; margin-top: 15px;",

        # 使用简单的 HTML 渲染（将 Markdown 转为 HTML）
        shiny::htmlOutput("ai_markdown_output")
      )
    )
  })

  # 导出按钮
  output$ai_export_buttons <- renderUI({
    result <- ai_result()

    if (is.null(result) || !result$success) {
      return(NULL)
    }

    # 检测 LaTeX 是否可用
    has_latex <- FALSE
    tryCatch({
      has_latex <- tinytex::tinytex_root() != "" || Sys.which("xelatex") != ""
    }, error = function(e) {
      has_latex <- Sys.which("xelatex") != ""
    })

    tags$div(
      style = "margin-top: 15px; padding: 15px; background: #f8f9fa; border-radius: 8px;",
      tags$h5("📥 导出报告"),
      tags$div(
        style = "display: flex; gap: 10px; flex-wrap: wrap;",
        downloadButton("download_ai_md", "📄 Markdown",
                       class = "btn-sm btn-outline-primary",
                       style = "border-radius: 5px;"),
        downloadButton("download_ai_html", "🌐 HTML",
                       class = "btn-sm btn-outline-success",
                       style = "border-radius: 5px;"),
        if (has_latex) {
          downloadButton("download_ai_pdf", "📕 PDF",
                         class = "btn-sm btn-outline-danger",
                         style = "border-radius: 5px;")
        } else {
          tags$button("📕 PDF (需安装LaTeX)",
                      class = "btn btn-sm btn-outline-secondary",
                      style = "border-radius: 5px; cursor: not-allowed;",
                      disabled = TRUE,
                      title = "请运行: install.packages('tinytex'); tinytex::install_tinytex()")
        }
      ),
      if (!has_latex) {
        tags$div(
          class = "alert alert-warning small",
          style = "margin-top: 10px; margin-bottom: 0; padding: 8px;",
          tags$strong("💡 提示: "),
          "PDF导出需要LaTeX。请运行以下命令安装：",
          tags$code("install.packages('tinytex'); tinytex::install_tinytex()"),
          tags$br(),
          "或直接使用 ",
          tags$strong("HTML导出"),
          "，然后在浏览器中 Ctrl+P 打印为PDF（推荐）"
        )
      } else {
        tags$p(class = "text-muted small", style = "margin-top: 10px; margin-bottom: 0;",
               "提示: HTML格式可在浏览器中打印为PDF，效果更佳")
      }
    )
  })

  # 渲染 Markdown 内容
  output$ai_markdown_output <- renderUI({
    result <- ai_result()
    if (is.null(result) || !result$success) return(NULL)

    # Markdown 到 HTML 转换
    content <- result$content

    # 分行处理
    lines <- strsplit(content, "\n")[[1]]
    html_lines <- character(length(lines))
    in_code_block <- FALSE
    in_list <- FALSE

    for (i in seq_along(lines)) {
      line <- lines[i]

      # 处理代码块
      if (grepl("^```", line)) {
        if (in_code_block) {
          html_lines[i] <- "</pre>"
          in_code_block <- FALSE
        } else {
          html_lines[i] <- "<pre style='background:#f0f0f0;padding:10px;overflow-x:auto;border-radius:4px;'>"
          in_code_block <- TRUE
        }
        next
      }

      # 代码块内保持原样
      if (in_code_block) {
        html_lines[i] <- line
        next
      }

      # 处理标题 (必须先检查更长的模式)
      if (grepl("^### ", line)) {
        if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
        html_lines[i] <- sub("^### (.*)", "<h4>\\1</h4>", line)
        next
      }
      if (grepl("^## ", line)) {
        if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
        html_lines[i] <- sub("^## (.*)", "<h3>\\1</h3>", line)
        next
      }
      if (grepl("^# ", line)) {
        if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
        html_lines[i] <- sub("^# (.*)", "<h2>\\1</h2>", line)
        next
      }

      # 处理列表项
      if (grepl("^- ", line)) {
        if (!in_list) {
          html_lines[i] <- paste0("<ul><li>", sub("^- ", "", line), "</li>")
          in_list <- TRUE
        } else {
          html_lines[i] <- paste0("<li>", sub("^- ", "", line), "</li>")
        }
        next
      }

      # 非列表行，如果之前在列表中则关闭
      if (in_list && !grepl("^- ", line) && nchar(trimws(line)) > 0) {
        html_lines[i] <- paste0("</ul>", line)
        in_list <- FALSE
        next
      }

      # 保持原行
      html_lines[i] <- line
    }

    # 关闭未闭合的标签
    if (in_code_block) html_lines <- c(html_lines, "</pre>")
    if (in_list) html_lines <- c(html_lines, "</ul>")

    # 重新组合
    content <- paste(html_lines, collapse = "\n")

    # 转换行内格式 (粗体、斜体、代码)
    content <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", content)
    content <- gsub("\\*(.+?)\\*", "<em>\\1</em>", content)
    content <- gsub("`([^`]+)`", "<code style='background:#f0f0f0;padding:2px 6px;border-radius:3px;'>\\1</code>", content)

    # 转换段落（连续两个换行）
    content <- gsub("\n\n+", "</p>\n<p>", content)

    # 转换单个换行
    content <- gsub("\n", "<br>\n", content)

    # 包装在段落中
    content <- paste0("<div style='line-height:1.8;'><p>", content, "</p></div>")

    shiny::HTML(content)
  })

  # 使用统计
  output$ai_usage_stats <- renderUI({
    history <- interpretation_history()

    if (nrow(history) == 0) {
      return(tags$p(class = "text-muted small", "暂无使用记录"))
    }

    total_tokens <- sum(history$tokens, na.rm = TRUE)

    tagList(
      tags$p(class = "small",
             tags$strong("📊 使用统计"),
             tags$br(),
             "总请求数: ", nrow(history),
             tags$br(),
             "总 Token: ", total_tokens
      )
    )
  })

  # 历史记录表格
  output$ai_history_table <- DT::renderDataTable({
    history <- interpretation_history()

    if (nrow(history) == 0) {
      df <- data.frame(
        时间 = character(),
        类型 = character(),
        模型 = character(),
        Token = integer(),
        stringsAsFactors = FALSE
      )
    } else {
      df <- data.frame(
        时间 = history$timestamp,
        类型 = history$type,
        模型 = history$model,
        Token = history$tokens,
        stringsAsFactors = FALSE
      )
    }

    DT::datatable(df,
                  options = list(
                    pageLength = 5,
                    language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Chinese.json')
                  ),
                  rownames = FALSE)
  })

  # 下载历史记录
  output$download_ai_history <- downloadHandler(
    filename = function() {
      paste0("ai_interpretation_history_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      history <- interpretation_history()
      write.csv(history, file, row.names = FALSE)
    }
  )

  # =====================================================
  # 导出 AI 解读结果
  # =====================================================

  # 下载 Markdown
  output$download_ai_md <- downloadHandler(
    filename = function() {
      paste0("AI解读报告_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".md")
    },
    content = function(file) {
      result <- ai_result()
      if (is.null(result) || !result$success) {
        showNotification("❌ 没有可导出的结果", type = "error")
        return()
      }

      # 生成带元数据的 Markdown
      header <- paste0(
        "---\n",
        "title: \"RNA-seq 差异表达分析 AI 解读报告\"\n",
        "date: \"", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\"\n",
        "model: \"", result$model, "\"\n",
        "tokens: ", result$tokens_used %||% "N/A", "\n",
        "generated_by: \"YuanSeq AI Interpretation\"\n",
        "---\n\n"
      )

      writeLines(paste0(header, result$content), file, useBytes = TRUE)
    }
  )

  # 下载 HTML（包含图片）
  output$download_ai_html <- downloadHandler(
    filename = function() {
      paste0("AI解读报告_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      result <- ai_result()
      if (is.null(result) || !result$success) {
        showNotification("❌ 没有可导出的结果", type = "error")
        return()
      }

      # 收集图片路径
      images <- list(
        volcano = isolate(plot_paths$volcano),
        heatmap = isolate(plot_paths$heatmap),
        kegg = isolate(plot_paths$kegg),
        go = isolate(plot_paths$go),
        tf = isolate(plot_paths$tf),
        pathway = isolate(plot_paths$pathway)
      )

      # 生成完整的 HTML 文档（包含图片）
      html_content <- markdown_to_html(result$content, result$model, result$tokens_used, images)
      writeLines(html_content, file, useBytes = TRUE)
    }
  )

  # 下载 PDF (需要 rmarkdown + pandoc)
  output$download_ai_pdf <- downloadHandler(
    filename = function() {
      paste0("AI解读报告_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      result <- ai_result()
      if (is.null(result) || !result$success) {
        showNotification("❌ 没有可导出的结果", type = "error")
        return()
      }

      # 检查 rmarkdown 是否可用
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        showNotification("❌ 导出PDF需要安装 rmarkdown 包: install.packages('rmarkdown')", type = "error", duration = 10)
        return()
      }

      # 创建临时 Rmd 文件
      temp_rmd <- tempfile(fileext = ".Rmd")

      rmd_content <- paste0(
        "---\n",
        "title: \"RNA-seq 差异表达分析 AI 解读报告\"\n",
        "date: \"", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\"\n",
        "output:\n",
        "  pdf_document:\n",
        "    toc: true\n",
        "    toc_depth: 3\n",
        "    number_sections: true\n",
        "    latex_engine: xelatex\n",
        "mainfont: \"SimSun\"\n",
        "---\n\n",
        "**模型**: ", result$model, "  \n",
        "**Token**: ", result$tokens_used %||% "N/A", "  \n",
        "**生成工具**: YuanSeq AI Interpretation\n\n",
        "---\n\n",
        result$content
      )

      writeLines(rmd_content, temp_rmd, useBytes = TRUE)

      tryCatch({
        # 渲染为 PDF
        rmarkdown::render(
          input = temp_rmd,
          output_file = file,
          output_format = rmarkdown::pdf_document(
            toc = TRUE,
            toc_depth = 3,
            number_sections = TRUE,
            latex_engine = "xelatex"
          ),
          quiet = TRUE
        )
      }, error = function(e) {
        showNotification(paste0("❌ PDF生成失败: ", e$message, "\n可能需要安装 LaTeX 或使用 HTML 导出"), type = "error", duration = 15)
      })
    }
  )

  # 返回结果（供其他模块使用）
  return(reactive({ ai_result() }))
}

# =====================================================
# 辅助函数：图片转 base64
# =====================================================

image_to_base64 <- function(file_path) {
  if (is.null(file_path) || !file.exists(file_path)) {
    return(NULL)
  }
  tryCatch({
    # 读取图片文件为二进制
    img_data <- readBin(file_path, "raw", file.info(file_path)$size)
    # 转换为base64
    base64_data <- base64enc::base64encode(file_path)
    # 根据文件扩展名确定MIME类型
    ext <- tolower(tools::file_ext(file_path))
    mime_type <- switch(ext,
                        png = "image/png",
                        jpg = "image/jpeg",
                        jpeg = "image/jpeg",
                        gif = "image/gif",
                        svg = "image/svg+xml",
                        "image/png")
    paste0("data:", mime_type, ";base64,", base64_data)
  }, error = function(e) {
    NULL
  })
}

# =====================================================
# 辅助函数：Markdown 转 HTML（支持图片）
# =====================================================

markdown_to_html <- function(content, model = "", tokens = NA, images = NULL) {
  # 分行处理
  lines <- strsplit(content, "\n")[[1]]
  html_lines <- character(length(lines))
  in_code_block <- FALSE
  in_list <- FALSE

  for (i in seq_along(lines)) {
    line <- lines[i]

    # 处理代码块
    if (grepl("^```", line)) {
      if (in_code_block) {
        html_lines[i] <- "</code></pre>"
        in_code_block <- FALSE
      } else {
        lang <- sub("^```(.*)$", "\\1", line)
        html_lines[i] <- paste0("<pre><code class=\"language-", lang, "\">")
        in_code_block <- TRUE
      }
      next
    }

    if (in_code_block) {
      html_lines[i] <- html_escape(line)
      next
    }

    # 处理标题
    if (grepl("^### ", line)) {
      if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
      html_lines[i] <- paste0("<h3>", sub("^### ", "", line), "</h3>")
      next
    }
    if (grepl("^## ", line)) {
      if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
      html_lines[i] <- paste0("<h2>", sub("^## ", "", line), "</h2>")
      next
    }
    if (grepl("^# ", line)) {
      if (in_list) { html_lines[i] <- "</ul>"; in_list <- FALSE }
      html_lines[i] <- paste0("<h1>", sub("^# ", "", line), "</h1>")
      next
    }

    # 处理列表
    if (grepl("^- ", line)) {
      if (!in_list) {
        html_lines[i] <- paste0("<ul><li>", sub("^- ", "", line), "</li>")
        in_list <- TRUE
      } else {
        html_lines[i] <- paste0("<li>", sub("^- ", "", line), "</li>")
      }
      next
    }

    if (in_list && !grepl("^- ", line) && nchar(trimws(line)) > 0) {
      html_lines[i] <- paste0("</ul>", line)
      in_list <- FALSE
      next
    }

    html_lines[i] <- line
  }

  if (in_code_block) html_lines <- c(html_lines, "</code></pre>")
  if (in_list) html_lines <- c(html_lines, "</ul>")

  content <- paste(html_lines, collapse = "\n")

  # 转换行内格式
  content <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", content)
  content <- gsub("\\*(.+?)\\*", "<em>\\1</em>", content)
  content <- gsub("`([^`]+)`", "<code>\\1</code>", content)

  # 转换段落
  content <- gsub("\n\n+", "</p>\n<p>", content)
  content <- gsub("\n", "<br>\n", content)

  # 生成图片HTML部分
  images_html <- ""
  if (!is.null(images)) {
    img_sections <- c()

    # 火山图
    if (!is.null(images$volcano) && file.exists(images$volcano)) {
      img_base64 <- image_to_base64(images$volcano)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>火山图 Volcano Plot</h3>\n',
          '  <img src="', img_base64, '" alt="火山图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    # 热图
    if (!is.null(images$heatmap) && file.exists(images$heatmap)) {
      img_base64 <- image_to_base64(images$heatmap)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>热图 Heatmap</h3>\n',
          '  <img src="', img_base64, '" alt="热图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    # KEGG通路图
    if (!is.null(images$kegg) && file.exists(images$kegg)) {
      img_base64 <- image_to_base64(images$kegg)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>KEGG 富集分析</h3>\n',
          '  <img src="', img_base64, '" alt="KEGG富集图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    # GO富集图
    if (!is.null(images$go) && file.exists(images$go)) {
      img_base64 <- image_to_base64(images$go)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>GO 富集分析</h3>\n',
          '  <img src="', img_base64, '" alt="GO富集图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    # TF活性图
    if (!is.null(images$tf) && file.exists(images$tf)) {
      img_base64 <- image_to_base64(images$tf)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>转录因子活性分析</h3>\n',
          '  <img src="', img_base64, '" alt="TF活性图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    # 通路活性图
    if (!is.null(images$pathway) && file.exists(images$pathway)) {
      img_base64 <- image_to_base64(images$pathway)
      if (!is.null(img_base64)) {
        img_sections <- c(img_sections, paste0(
          '<div class="figure">\n',
          '  <h3>通路活性分析</h3>\n',
          '  <img src="', img_base64, '" alt="通路活性图" style="max-width: 100%; height: auto;">\n',
          '</div>\n'
        ))
      }
    }

    if (length(img_sections) > 0) {
      images_html <- paste0(
        '\n<hr>\n',
        '<h2>分析图表</h2>\n',
        paste(img_sections, collapse = "\n")
      )
    }
  }

  # 生成完整 HTML 文档
  html_template <- paste0(
    '<!DOCTYPE html>\n',
    '<html lang="zh-CN">\n',
    '<head>\n',
    '  <meta charset="UTF-8">\n',
    '  <meta name="viewport" content="width=device-width, initial-scale=1.0">\n',
    '  <title>RNA-seq 差异表达分析 AI 解读报告</title>\n',
    '  <style>\n',
    '    body { font-family: "Microsoft YaHei", "SimSun", Arial, sans-serif; line-height: 1.8; max-width: 900px; margin: 0 auto; padding: 20px; color: #333; }\n',
    '    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }\n',
    '    h2 { color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; margin-top: 30px; }\n',
    '    h3 { color: #7f8c8d; margin-top: 20px; }\n',
    '    pre { background: #f8f9fa; padding: 15px; border-radius: 5px; overflow-x: auto; border-left: 4px solid #3498db; }\n',
    '    code { background: #f8f9fa; padding: 2px 6px; border-radius: 3px; font-family: Consolas, monospace; }\n',
    '    pre code { background: none; padding: 0; }\n',
    '    ul { padding-left: 20px; }\n',
    '    li { margin: 8px 0; }\n',
    '    .meta { background: #ecf0f1; padding: 15px; border-radius: 5px; margin-bottom: 20px; }\n',
    '    .meta span { margin-right: 20px; }\n',
    '    .figure { margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 8px; text-align: center; }\n',
    '    .figure img { box-shadow: 0 2px 8px rgba(0,0,0,0.1); border-radius: 4px; }\n',
    '    @media print { body { max-width: 100%; } .figure { page-break-inside: avoid; } }\n',
    '  </style>\n',
    '</head>\n',
    '<body>\n',
    '  <h1>RNA-seq 差异表达分析 AI 解读报告</h1>\n',
    '  <div class="meta">\n',
    '    <span><strong>模型:</strong> ', model, '</span>\n',
    '    <span><strong>Token:</strong> ', tokens %||% 'N/A', '</span>\n',
    '    <span><strong>日期:</strong> ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</span>\n',
    '    <span><strong>工具:</strong> YuanSeq</span>\n',
    '  </div>\n',
    '  <hr>\n',
    '  <div class="content">\n',
    '  <p>', content, '</p>\n',
    '  </div>\n',
    images_html, '\n',
    '  <hr>\n',
    '  <p style="text-align: center; color: #7f8c8d; font-size: 12px;">Generated by YuanSeq AI Interpretation</p>\n',
    '</body>\n',
    '</html>'
  )

  return(html_template)
}

# HTML 转义
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x)
  x <- gsub("<", "&lt;", x)
  x <- gsub(">", "&gt;", x)
  x <- gsub("\"", "&quot;", x)
  x <- gsub("'", "&#39;", x)
  return(x)
}
