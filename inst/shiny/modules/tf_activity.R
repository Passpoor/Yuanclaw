# =====================================================
# 转录因子活性模块
# =====================================================

tf_activity_server <- function(input, output, session, deg_results) {

  # 缓存 CollecTRI 网络
  collectri_net <- reactive({
    req(input$data_source)

    if (input$data_source == "counts") {
      species_code <- input$species_select
    } else {
      species_code <- input$deg_species
    }

    organism_code <- if(species_code == "Hs") 'human' else 'mouse'
    current_file <- paste0("collectri_", organism_code, ".rds")

    if (file.exists(current_file)) {
      showNotification(paste0("正在从本地加载 CollecTRI (", organism_code, ")..."), type = "message", duration = 3)
      return(readRDS(current_file))
    }

    showNotification(paste0("本地文件不存在，正在从网上下载 CollecTRI (", organism_code, ")..."), type = "warning", duration = 5)

    tryCatch({
      net <- decoupleR::get_collectri(organism = organism_code, split_complexes = FALSE)
      saveRDS(net, current_file)
      showNotification(paste0("CollecTRI (", organism_code, ") 下载并保存成功!"), type = "message", duration = 3)
      return(net)
    }, error = function(e) {
      showNotification(paste("下载 CollecTRI 网络失败:", e$message), type = "error")
      return(NULL)
    })
  })

  # 运行 TF 活性分析 - 支持多算法
  tf_activity_results <- eventReactive(input$run_tf_activity, {
    req(deg_results(), collectri_net())

    # 🆕 获取选择的算法
    method <- input$tf_method

    showNotification(paste("正在推断转录因子活性 (算法:", toupper(method), ")..."), type = "message")

    # 从deg_results中提取差异分析结果
    deg_data <- deg_results()
    deg_res <- deg_data$deg_df
    net <- collectri_net()

    # 1. 构造输入矩阵 (t-like 统计量)
    stats_df <- deg_res %>%
      filter(!is.na(SYMBOL), !is.na(t_stat)) %>%  # 🌟 改为SYMBOL

      group_by(SYMBOL) %>%
      filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
      ungroup() %>%
      distinct(SYMBOL, .keep_all = TRUE)

    if (nrow(stats_df) < 5) {
      showNotification(
        paste0("TF 分析失败: 用于 TF 推断的基因数量 (", nrow(stats_df), ") 不足，请检查数据和阈值。"),
        type = "error",
        duration = 15
      )
      return(NULL)
    }

    # 2. ID 兼容性检查
    input_genes <- stats_df$SYMBOL  # 🌟 改为SYMBOL
    net_targets <- unique(net$target)
    shared_genes <- intersect(input_genes, net_targets)

    if(length(shared_genes) < input$tf_min_size) {
      showNotification(
        paste0("TF 分析失败: 共享靶基因数量 (", length(shared_genes), ") 小于最小要求 (", input$tf_min_size, ")。"),
        type = "error",
        duration = 15
      )
      return(NULL)
    }

    stats_df_filtered <- stats_df %>% filter(SYMBOL %in% shared_genes)  # 🌟 改为SYMBOL

    # 🔥 关键修复：移除NA和Inf值
    stats_df_clean <- stats_df_filtered %>%
      filter(!is.na(t_stat)) %>%           # 移除NA
      filter(is.finite(t_stat)) %>%        # 移除Inf和-Inf
      filter(t_stat != 0)                   # 移除0值（可选）

    cat(sprintf("📊 TF分析: 原始 %d 基因 -> 清洗后 %d 基因\n",
                nrow(stats_df_filtered), nrow(stats_df_clean)))

    if (nrow(stats_df_clean) < 5) {
      showNotification(
        paste0("TF 分析失败: 清洗后的有效基因数量 (", nrow(stats_df_clean), ") 不足，请检查数据质量。"),
        type = "error",
        duration = 15
      )
      return(NULL)
    }

    mat_input <- stats_df_clean %>%
      select(SYMBOL, t_stat) %>%  # 🌟 改为SYMBOL
      column_to_rownames(var = "SYMBOL") %>%  # 🌟 改为SYMBOL
      as.matrix()

    # 检查矩阵是否还有NA或Inf
    if (any(is.na(mat_input)) || any(!is.finite(mat_input))) {
      cat("⚠️ 警告: 矩阵中仍有NA或Inf值\n")
      mat_input <- mat_input[is.finite(rowSums(mat_input)), ]
      mat_input <- mat_input[, is.finite(colSums(mat_input))]
    }

    # 3. 🆕 根据选择的算法运行
    tryCatch({
      contrast_acts <- switch(method,
        "ulm" = {
          decoupleR::run_ulm(
            mat = mat_input,
            net = net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$tf_min_size
          )
        },
        "mlm" = {
          decoupleR::run_mlm(
            mat = mat_input,
            net = net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$tf_min_size
          )
        },
        "wmean" = {
          decoupleR::run_wmean(
            mat = mat_input,
            net = net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$tf_min_size
          )
        },
        "wsum" = {
          decoupleR::run_wsum(
            mat = mat_input,
            net = net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$tf_min_size
          )
        },
        {
          # 默认使用ULM
          decoupleR::run_ulm(
            mat = mat_input,
            net = net,
            .source = 'source',
            .target = 'target',
            .mor = 'mor',
            minsize = input$tf_min_size
          )
        }
      )

      # 确保返回的是数据框
      if (is.data.frame(contrast_acts)) {
        result_df <- contrast_acts
      } else if (is.list(contrast_acts)) {
        # 某些算法可能返回列表
        if ("statistic" %in% names(contrast_acts)) {
          result_df <- contrast_acts$statistic
        } else {
          result_df <- as.data.frame(contrast_acts)
        }
      } else {
        result_df <- as.data.frame(contrast_acts)
      }

      # 4. 添加排名信息
      result_df <- result_df %>%
        mutate(rnk = NA)

      msk_pos <- result_df$score > 0
      result_df[msk_pos, 'rnk'] <- rank(-result_df[msk_pos, 'score'])

      msk_neg <- result_df$score < 0
      result_df[msk_neg, 'rnk'] <- rank(-abs(result_df[msk_neg, 'score']))

      # 添加方法信息
      result_df$method <- toupper(method)

      showNotification(paste("转录因子活性推断完成! (算法:", toupper(method), ")"), type = "message")

      return(result_df)

    }, error = function(e) {
      error_msg <- e$message

      # 🔥 为MLM共线性错误提供特殊提示
      if (grepl("colinear", error_msg, ignore.case = TRUE)) {
        showNotification(
          paste("TF 活性分析失败 (", toupper(method), "): 检测到共线性问题。",
                "\n建议: 请尝试使用ULM、WMEAN或WSUM算法代替MLM算法。",
                "\nMLM算法对数据质量要求较高，容易出现共线性问题。"),
          type = "error",
          duration = 10
        )
      } else {
        showNotification(paste("TF 活性分析失败 (", toupper(method), "):", error_msg), type = "error")
      }

      return(NULL)
    })
  })

  # 绘制 TF 活性柱状图
  output$tf_activity_bar_plot <- renderPlot({
    req(tf_activity_results())

    df_acts <- tf_activity_results()

    n_tfs <- input$tf_top_n

    tfs_to_plot <- df_acts %>%
      arrange(rnk) %>%
      head(n_tfs) %>%
      pull(source)

    f_contrast_acts <- df_acts %>%
      filter(source %in% tfs_to_plot)

    txt_col <- if(input$theme_toggle) "white" else "black"
    grid_col <- if(input$theme_toggle) "#444444" else "#cccccc"

    ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) +
      geom_bar(aes(fill = score), stat = "identity") +
      scale_fill_gradient2(
        low = input$tf_inactive_col,
        high = input$tf_active_col,
        mid = "whitesmoke",
        midpoint = 0,
        name = "TF 活性分数"
      ) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = txt_col) +
      theme_minimal() +
      labs(x = "转录因子 (TFs)", y = "活性分数 (ULM Score)",
           title = paste("Top", n_tfs, "转录因子活性变化 (T vs C)")) +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        axis.title = element_text(color = txt_col, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = txt_col),
        axis.text.y = element_text(size = 10, face = "bold", color = txt_col),
        legend.text = element_text(color = txt_col),
        legend.title = element_text(color = txt_col),
        axis.line = element_line(color = txt_col),
        panel.grid.major = element_line(color = grid_col),
        panel.grid.minor = element_line(color = grid_col)
      )
  })

  output$download_tf_results <- downloadHandler(
    filename = function() {
      paste0("TF_Activity_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(tf_activity_results())
      df <- tf_activity_results() %>%
        select(source, score, p_value, rnk) %>%
        rename(TF = source, Score = score, P.Value = p_value, Rank = rnk)
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$tf_activity_table <- DT::renderDataTable({
    req(tf_activity_results())

    df <- tf_activity_results() %>%
      select(source, score, p_value, rnk) %>%
      rename(TF = source, Score = score, P.Value = p_value, Rank = rnk) %>%
      arrange(Rank)

    DT::datatable(df, selection = 'single', options = list(scrollX=T, pageLength=10), rownames=F) %>%
      formatRound(c("Score", "P.Value"), 4)
  })

  # =====================================================
  # 修复：selected_tf_targets 函数 - 修复select()错误
  # =====================================================
  selected_tf_targets <- reactive({
    req(tf_activity_results(), collectri_net(), deg_results())

    selected_row <- input$tf_activity_table_rows_selected
    if (length(selected_row) == 0) {
      return(NULL)
    }

    # 🔥 关键修复：添加类型检查
    tf_res <- tf_activity_results()
    if (!is.data.frame(tf_res)) {
      showNotification("TF结果格式错误: 不是数据框", type = "error")
      return(NULL)
    }

    net <- collectri_net()
    if (!is.data.frame(net)) {
      showNotification("CollecTRI网络格式错误: 不是数据框", type = "error")
      return(NULL)
    }

    deg_res <- deg_results()
    if (!is.list(deg_res) || is.null(deg_res$deg_df)) {
      showNotification("差异分析结果格式错误", type = "error")
      return(NULL)
    }

    tf_name <- tf_res %>%
      arrange(rnk) %>%
      slice(selected_row) %>%
      pull(source)

    tf_net <- net %>%
      filter(source == tf_name) %>%
      select(target, mor) %>%
      rename(SYMBOL = target, Mode_of_Regulation = mor)

    # 🔥 关键修复：使用 $deg_df 提取数据框
    deg_data <- deg_res$deg_df %>%
      select(SYMBOL, log2FoldChange, pvalue, padj, t_stat, Status)

    final_table <- tf_net %>%
      left_join(deg_data, by = "SYMBOL") %>%  # 🌟 改为SYMBOL
      mutate(
        t_stat_clean = tidyr::replace_na(t_stat, 0),
        mor_clean = tidyr::replace_na(Mode_of_Regulation, 0)
      ) %>%
      mutate(Is_DE = ifelse(Status != "Not DE", "Yes", "No") ) %>%
      # 🔥 关键：添加 Match_Status 列
      mutate(
        Predicted_Change = case_when(
          Mode_of_Regulation > 0 ~ "Activator (Up)",
          Mode_of_Regulation < 0 ~ "Repressor (Down)",
          TRUE ~ "Unknown"
        ),
        Actual_Change = case_when(
          log2FoldChange > 0 ~ "Up",
          log2FoldChange < 0 ~ "Down",
          TRUE ~ "Not DE"
        ),
        Match_Status = case_when(
          Predicted_Change == "Activator (Up)" & Actual_Change == "Up" ~ "Consistent",
          Predicted_Change == "Repressor (Down)" & Actual_Change == "Down" ~ "Consistent",
          Predicted_Change == "Activator (Up)" & Actual_Change == "Down" ~ "Inconsistent",
          Predicted_Change == "Repressor (Down)" & Actual_Change == "Up" ~ "Inconsistent",
          TRUE ~ "Neutral/Unknown"
        )
      ) %>%
      select(SYMBOL, Mode_of_Regulation, Is_DE, Status, log2FoldChange, t_stat, pvalue, padj,
             Predicted_Change, Actual_Change, Match_Status) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      # 移除临时列
      select(-Predicted_Change, -Actual_Change)

    return(list(tf_name = tf_name, data = final_table))
  })  # 🌟 这里添加了缺失的右括号

  output$tf_target_table <- DT::renderDataTable({
    req(selected_tf_targets())

    data_list <- selected_tf_targets()

    DT::datatable(
      data_list$data,
      caption = paste0("TF: ", data_list$tf_name, " 的靶基因"),
      options = list(scrollX=T, pageLength=10), rownames=F
    ) %>%
      formatRound(c("log2FoldChange", "t_stat", "pvalue", "padj"), 4)
  })

  # 🆕 一致性统计
  output$tf_consistency_summary <- renderText({
    req(selected_tf_targets())

    df <- selected_tf_targets()$data

    # 计算一致性统计
    n_consistent <- sum(df$Match_Status == "Consistent", na.rm = TRUE)
    n_inconsistent <- sum(df$Match_Status == "Inconsistent", na.rm = TRUE)
    n_neutral <- sum(df$Match_Status == "Neutral/Unknown", na.rm = TRUE)
    n_total <- n_consistent + n_inconsistent + n_neutral

    if (n_total > 0) {
      pct_consistent <- round(100 * n_consistent / n_total, 1)
      pct_inconsistent <- round(100 * n_inconsistent / n_total, 1)
      pct_neutral <- round(100 * n_neutral / n_total, 1)

      paste0(
        "📊 调控一致性统计 | ",
        "✅ 一致: ", n_consistent, " (", pct_consistent, "%) | ",
        "❌ 不一致: ", n_inconsistent, " (", pct_inconsistent, "%) | ",
        "⚪ 未知: ", n_neutral, " (", pct_neutral, "%)"
      )
    } else {
      "📊 无数据"
    }
  })

  output$tf_target_plot <- renderPlot({
    req(selected_tf_targets())

    # 🔍 获取自定义参数
    point_size <- input$tf_scatter_point_size
    if (is.null(point_size)) point_size <- 3

    point_alpha <- input$tf_scatter_alpha
    if (is.null(point_alpha)) point_alpha <- 0.7

    label_size <- input$tf_scatter_label_size
    if (is.null(label_size)) label_size <- 3

    n_labels <- input$tf_scatter_n_labels
    if (is.null(n_labels)) n_labels <- 15

    data_list <- selected_tf_targets()
    df <- data_list$data
    tf_name <- data_list$tf_name

    # 🔥 Match_Status 已经在 selected_tf_targets() 中计算好了
    # 不需要重复计算

    txt_col <- if(input$theme_toggle) "white" else "black"
    grid_col <- if(input$theme_toggle) "#444444" else "#cccccc"

    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = Match_Status), size = point_size, alpha = point_alpha) +
      scale_color_manual(
        values = c("Consistent" = "#2ecc71", "Inconsistent" = "#e74c3c", "Neutral/Unknown" = "#95a5a6"),
        name = "调控一致性"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = txt_col, alpha = 0.7) +
      geom_hline(yintercept = -log10(input$pval_cutoff), linetype = "dotted", color = txt_col, alpha = 0.7) +
      labs(
        title = paste("TF:", tf_name, "的靶基因差异表达 (Target Genes DE Plot)"),
        x = "log2(Fold Change)",
        y = "-log10(P Value)"
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        axis.title = element_text(color = txt_col, face = "bold"),
        axis.text = element_text(color = txt_col),
        legend.text = element_text(color = txt_col),
        legend.title = element_text(color = txt_col),
        axis.line = element_line(color = txt_col),
        panel.grid.major = element_line(color = grid_col),
        panel.grid.minor = element_line(color = grid_col)
      )

    # 🆕 添加top基因的标签 - 使用用户自定义数量
    if (n_labels > 0) {
      top_genes <- df %>%
        filter(!is.na(pvalue)) %>%
        arrange(pvalue) %>%
        head(min(n_labels, nrow(df)))  # 使用用户自定义的数量

      if (nrow(top_genes) > 0) {
        # 为标签添加偏移，避免重叠
        top_genes <- top_genes %>%
          mutate(
            label_x = log2FoldChange + ifelse(log2FoldChange > 0, 0.2, -0.2),
            label_y = -log10(pvalue) + 0.5
          )

        p <- p +
          geom_text(
            data = top_genes,
            aes(x = label_x, y = label_y, label = SYMBOL),
            size = label_size,  # 使用用户自定义的标签大小
            color = txt_col,
            fontface = "bold",
            check_overlap = TRUE,
            vjust = 0.5,
            hjust = ifelse(top_genes$log2FoldChange > 0, 0, 1)
          )
      }
    }

    print(p)
  })

  # 🆕 交互式TF靶基因网络可视化 - 支持拖动
  output$tf_network_plot_interactive <- renderPlotly({
    req(selected_tf_targets())

    # 🔥 强制响应所有自定义选项的变化 - 直接引用input
    input$tf_network_node_size
    input$tf_network_label_size
    input$tf_tf_node_col
    input$tf_consistent_act_col
    input$tf_consistent_rep_col
    input$tf_inconsistent_act_col
    input$tf_inconsistent_rep_col
    input$tf_neutral_col
    input$theme_toggle

    # 🔍 确保input值存在
    node_size_mult <- input$tf_network_node_size
    if (is.null(node_size_mult)) node_size_mult <- 1

    # 🔥 获取标签大小设置
    label_size <- input$tf_network_label_size
    if (is.null(label_size)) label_size <- 3.5

    data_list <- selected_tf_targets()
    df <- data_list$data
    tf_name <- data_list$tf_name

    # 选择top靶基因（按pvalue和log2FC）
    top_targets <- df %>%
      filter(!is.na(pvalue)) %>%
      arrange(pvalue, abs(log2FoldChange)) %>%
      head(min(30, nrow(df))) %>%
      mutate(
        node_type = "target",
        # 根据调控模式设置边颜色
        edge_color = ifelse(Mode_of_Regulation > 0, "#e74c3c", "#3498db"),
        # 根据DE状态设置节点大小
        node_size = ifelse(Is_DE == "Yes", 8, 5)
      )

    if (nrow(top_targets) < 2) {
      return(NULL)
    }

    # 靶基因节点环绕TF排列
    n_targets <- nrow(top_targets)
    angles <- seq(0, 2*pi, length.out = n_targets + 1)[1:n_targets]
    top_targets$x <- 2 * cos(angles)
    top_targets$y <- 2 * sin(angles)

    # 🆕 根据用户自定义颜色设置节点颜色
    top_targets$color <- ifelse(top_targets$Match_Status == "Consistent",
                                 ifelse(top_targets$Mode_of_Regulation > 0,
                                        input$tf_consistent_act_col,
                                        input$tf_consistent_rep_col),
                                 ifelse(top_targets$Match_Status == "Inconsistent",
                                        ifelse(top_targets$Mode_of_Regulation > 0,
                                               input$tf_inconsistent_act_col,
                                               input$tf_inconsistent_rep_col),
                                        input$tf_neutral_col))

    # 创建边的数据（每条边单独一行）
    edges_list <- list()
    for (i in 1:n_targets) {
      edges_list[[i]] <- list(
        x = c(0, top_targets$x[i], NA),  # TF -> 靶基因
        y = c(0, top_targets$y[i], NA),
        color = top_targets$edge_color[i],
        linetype = ifelse(top_targets$Match_Status[i] == "Consistent", "solid", "dashed")
      )
    }

    # 🔥 修复：简化代码，直接使用edges_list，不需要创建中间traces

    # 创建TF节点数据
    tf_node_data <- data.frame(x = 0, y = 0)

    # 使用subplot的简单模式 - share=TRUE共享坐标轴
    p <- plot_ly() %>%
      # 添加TF节点
      add_trace(
        data = tf_node_data,
        x = ~x, y = ~y,
        type = 'scatter',
        mode = 'markers+text',  # 🆕 添加text模式
        name = tf_name,
        text = ~tf_name,
        textfont = list(
          size = label_size,
          color = if(input$theme_toggle) "white" else "black"
        ),
        textposition = 'top center',
        marker = list(
          size = 12 * node_size_mult,
          color = input$tf_tf_node_col,
          line = list(color = 'white', width = 2)
        ),
        hoverinfo = 'text',
        showlegend = FALSE
      ) %>%
      # 添加靶基因节点
      add_trace(
        data = top_targets,
        x = ~x, y = ~y,
        type = 'scatter',
        mode = 'markers+text',  # 🆕 添加text模式
        name = 'Target Genes',
        text = ~SYMBOL,  # 🆕 显示基因名称
        textfont = list(
          size = label_size * 0.8,  # 稍小一点
          color = if(input$theme_toggle) "white" else "black"
        ),
        textposition = 'top center',
        hovertext = ~paste(
          "Gene:", SYMBOL, "<br>",
          "log2FC:", round(log2FoldChange, 3), "<br>",
          "p-value:", format(pvalue, scientific = TRUE, digits = 3), "<br>",
          "Status:", Match_Status
        ),
        marker = list(
          size = ~node_size * node_size_mult,
          color = ~color,
          line = list(color = 'white', width = 1)
        ),
        hoverinfo = 'text',
        showlegend = FALSE
      ) %>%
      # 添加所有边
      layout(
        title = paste("TF:", tf_name, "的靶基因调控网络（可拖动节点）"),
        showlegend = FALSE,
        xaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = c(-3.5, 3.5)
        ),
        yaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          scaleanchor = "x",
          range = c(-3.5, 3.5)
        ),
        plot_bgcolor = if(input$theme_toggle) "#2b2b2b" else "white",
        paper_bgcolor = if(input$theme_toggle) "#1a1a1a" else "white",
        font = list(color = if(input$theme_toggle) "white" else "black"),
        hovermode = 'closest',
        dragmode = 'move'  # 🔥 启用拖动模式
      )

    # 逐个添加边（直接从edges_list添加，避免数据提取问题）
    for (i in 1:length(edges_list)) {
      edge_data <- data.frame(
        x = edges_list[[i]]$x,
        y = edges_list[[i]]$y
      )
      p <- p %>% add_trace(
        data = edge_data,
        x = ~x, y = ~y,
        type = 'scatter',
        mode = 'lines',
        line = list(color = edges_list[[i]]$color, width = 2),
        showlegend = FALSE,
        hoverinfo = 'skip',
        inherit = FALSE
      )
    }

    p
  })

  # 🆕 TF靶基因网络可视化（原有静态图）
  output$tf_network_plot <- renderPlot({
    req(selected_tf_targets())

    # 🔍 确保input值存在
    node_size_mult <- input$tf_network_node_size
    if (is.null(node_size_mult)) node_size_mult <- 1

    label_size <- input$tf_network_label_size
    if (is.null(label_size)) label_size <- 3.5

    data_list <- selected_tf_targets()
    df <- data_list$data
    tf_name <- data_list$tf_name

    # 选择top靶基因（按pvalue和log2FC）
    top_targets <- df %>%
      filter(!is.na(pvalue)) %>%
      arrange(pvalue, abs(log2FoldChange)) %>%
      head(min(30, nrow(df))) %>%
      mutate(
        node_type = "target",
        # 根据调控模式设置颜色
        edge_color = ifelse(Mode_of_Regulation > 0, "#e74c3c", "#3498db"),  # 激活红色，抑制蓝色
        # 根据DE状态设置节点大小
        node_size = ifelse(Is_DE == "Yes", 8, 5)
      )

    if (nrow(top_targets) < 2) {
      # 数据太少，返回空图
      plot.new()
      title(main = paste("TF:", tf_name, "- 数据不足以绘制网络图"))
      return()
    }

    # 创建网络节点数据
    # TF节点在中心
    tf_node <- data.frame(
      name = tf_name,
      node_type = "tf",
      x = 0,
      y = 0,
      node_size = 12,
      color = input$tf_tf_node_col  # 🆕 使用自定义颜色
    )

    # 靶基因节点环绕TF排列
    n_targets <- nrow(top_targets)
    angles <- seq(0, 2*pi, length.out = n_targets + 1)[1:n_targets]

    top_targets$x <- 2 * cos(angles)
    top_targets$y <- 2 * sin(angles)

    # 🆕 根据用户自定义颜色设置节点颜色
    top_targets$color <- ifelse(top_targets$Match_Status == "Consistent",
                                 ifelse(top_targets$Mode_of_Regulation > 0,
                                        input$tf_consistent_act_col,  # 一致-激活
                                        input$tf_consistent_rep_col), # 一致-抑制
                                 ifelse(top_targets$Match_Status == "Inconsistent",
                                        ifelse(top_targets$Mode_of_Regulation > 0,
                                               input$tf_inconsistent_act_col,  # 不一致-激活
                                               input$tf_inconsistent_rep_col), # 不一致-抑制
                                        input$tf_neutral_col))  # 未知

    # 合并节点
    all_nodes <- rbind(
      tf_node[, c("name", "x", "y", "node_size", "color")],
      top_targets[, c("SYMBOL", "x", "y", "node_size", "color")] %>%
        rename(name = SYMBOL)
    )

    # 创建边数据 - 🔥 修复：添加坐标列
    edges <- data.frame(
      x = rep(0, n_targets),           # TF的x坐标（中心）
      y = rep(0, n_targets),           # TF的y坐标（中心）
      xend = top_targets$x,            # 靶基因的x坐标
      yend = top_targets$y,            # 靶基因的y坐标
      color = top_targets$edge_color,
      linetype = ifelse(top_targets$Match_Status == "Consistent", "solid", "dashed")
    )

    txt_col <- if(input$theme_toggle) "white" else "black"

    # 🔍 测试：强制使用固定值来测试
    # node_size_mult <- 2  # 测试：强制使用2倍
    # label_size <- 5       # 测试：强制使用5
    # cat("🔍🔍🔍 测试模式：强制使用节点大小=2，标签大小=5\n")

    # 🆕 应用节点大小倍数
    # 🔍 调试输出
    cat(sprintf("🔍 调试: tf_network_node_size = %s (使用值: %s)\n",
                input$tf_network_node_size, node_size_mult))
    cat(sprintf("🔍 调试: tf_network_label_size = %s (使用值: %s)\n",
                input$tf_network_label_size, label_size))

    all_nodes$node_size_scaled <- all_nodes$node_size * node_size_mult

    cat(sprintf("🔍 调试: 节点大小缩放: %s -> %s\n",
                all_nodes$node_size[1], all_nodes$node_size_scaled[1]))

    # 绘制网络图
    p <- ggplot() +
      # 绘制边
      geom_segment(
        data = edges,
        aes(x = x, xend = xend, y = y, yend = yend, color = color, linetype = linetype),
        linewidth = 0.8,  # 🆕 使用linewidth代替size (ggplot2 3.4.0+)
        alpha = 0.6
      ) +
      # 绘制节点
      geom_point(
        data = all_nodes,
        aes(x = x, y = y, size = node_size_scaled, color = color),
        alpha = 0.9
      ) +
      # 添加基因标签
      geom_text(
        data = all_nodes,
        aes(x = x, y = y, label = name),
        size = label_size,  # 🆕 使用变量
        fontface = "bold",
        color = txt_col,
        vjust = ifelse(all_nodes$y > 0, -0.5, 1.5),
        check_overlap = TRUE
      ) +
      scale_color_identity() +
      scale_linetype_identity() +
      scale_size_identity() +
      labs(
        title = paste("TF:", tf_name, "的靶基因调控网络"),
        subtitle = paste("Top", n_targets, "靶基因 | 红色=激活, 蓝色=抑制, 实线=一致, 虚线=不一致")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = txt_col, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
      ) +
      coord_equal() +
      xlim(-3, 3) +
      ylim(-3, 3)

    print(p)
  })

  # 🆕 SVG导出功能
  output$download_tf_network_svg <- downloadHandler(
    filename = function() {
      req(selected_tf_targets())
      tf_name <- selected_tf_targets()$tf_name
      paste0("TF_Network_", tf_name, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      req(selected_tf_targets())

      # 重新生成网络图（与上面相同）
      data_list <- selected_tf_targets()
      df <- data_list$data
      tf_name <- data_list$tf_name

      top_targets <- df %>%
        filter(!is.na(pvalue)) %>%
        arrange(pvalue, abs(log2FoldChange)) %>%
        head(min(30, nrow(df))) %>%
        mutate(
          node_type = "target",
          edge_color = ifelse(Mode_of_Regulation > 0, "#e74c3c", "#3498db"),
          node_size = ifelse(Is_DE == "Yes", 8, 5)
        )

      if (nrow(top_targets) < 2) return()

      tf_node <- data.frame(
        name = tf_name, node_type = "tf", x = 0, y = 0,
        node_size = 12, color = input$tf_tf_node_col
      )

      n_targets <- nrow(top_targets)
      angles <- seq(0, 2*pi, length.out = n_targets + 1)[1:n_targets]
      top_targets$x <- 2 * cos(angles)
      top_targets$y <- 2 * sin(angles)
      top_targets$color <- ifelse(top_targets$Match_Status == "Consistent",
                                   ifelse(top_targets$Mode_of_Regulation > 0,
                                          input$tf_consistent_act_col,
                                          input$tf_consistent_rep_col),
                                   ifelse(top_targets$Match_Status == "Inconsistent",
                                          ifelse(top_targets$Mode_of_Regulation > 0,
                                                 input$tf_inconsistent_act_col,
                                                 input$tf_inconsistent_rep_col),
                                          input$tf_neutral_col))

      all_nodes <- rbind(
        tf_node[, c("name", "x", "y", "node_size", "color")],
        top_targets[, c("SYMBOL", "x", "y", "node_size", "color")] %>% rename(name = SYMBOL)
      )

      edges <- data.frame(
        x = rep(0, n_targets), y = rep(0, n_targets),
        xend = top_targets$x, yend = top_targets$y,
        color = top_targets$edge_color,
        linetype = ifelse(top_targets$Match_Status == "Consistent", "solid", "dashed")
      )

      all_nodes$node_size_scaled <- all_nodes$node_size * input$tf_network_node_size

      txt_col <- if(input$theme_toggle) "white" else "black"

      p <- ggplot() +
        geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend,
                                        color = color, linetype = linetype),
                    linewidth = 0.8, alpha = 0.6) +
        geom_point(data = all_nodes, aes(x = x, y = y, size = node_size_scaled, color = color), alpha = 0.9) +
        geom_text(data = all_nodes, aes(x = x, y = y, label = name),
                  size = input$tf_network_label_size, fontface = "bold",
                  color = txt_col, vjust = ifelse(all_nodes$y > 0, -0.5, 1.5), check_overlap = TRUE) +
        scale_color_identity() + scale_linetype_identity() + scale_size_identity() +
        labs(title = paste("TF:", tf_name, "的靶基因调控网络"),
             subtitle = paste("Top", n_targets, "靶基因")) +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white", colour = NA),
              plot.title = element_text(face = "bold", hjust = 0.5),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), panel.grid = element_blank()) +
        coord_equal() + xlim(-3, 3) + ylim(-3, 3)

      # 保存为SVG
      svg(file, width = 10, height = 8)
      print(p)
      dev.off()
    },
    contentType = "image/svg+xml"
  )

  output$download_tf_scatter_svg <- downloadHandler(
    filename = function() {
      req(selected_tf_targets())
      tf_name <- selected_tf_targets()$tf_name
      paste0("TF_Scatter_", tf_name, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      req(selected_tf_targets())

      # 重新生成散点图 - 使用用户自定义设置
      data_list <- selected_tf_targets()
      df <- data_list$data
      tf_name <- data_list$tf_name

      # 🔍 获取自定义参数
      point_size <- input$tf_scatter_point_size
      if (is.null(point_size)) point_size <- 3

      point_alpha <- input$tf_scatter_alpha
      if (is.null(point_alpha)) point_alpha <- 0.7

      label_size <- input$tf_scatter_label_size
      if (is.null(label_size)) label_size <- 3

      n_labels <- input$tf_scatter_n_labels
      if (is.null(n_labels)) n_labels <- 15

      txt_col <- "black"  # SVG使用黑色文本
      grid_col <- "#cccccc"

      p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = Match_Status), size = point_size, alpha = point_alpha) +
        scale_color_manual(
          values = c("Consistent" = "#2ecc71", "Inconsistent" = "#e74c3c", "Neutral/Unknown" = "#95a5a6"),
          name = "调控一致性"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = txt_col, alpha = 0.7) +
        geom_hline(yintercept = -log10(input$pval_cutoff), linetype = "dotted", color = txt_col, alpha = 0.7) +
        labs(title = paste("TF:", tf_name, "的靶基因差异表达"),
             x = "log2(Fold Change)", y = "-log10(P Value)") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold")
        )

      # 添加基因标签 - 使用用户自定义数量
      if (n_labels > 0) {
        top_genes <- df %>%
          filter(!is.na(pvalue)) %>%
          arrange(pvalue) %>%
          head(min(n_labels, nrow(df)))

        if (nrow(top_genes) > 0) {
          top_genes <- top_genes %>%
            mutate(
              label_x = log2FoldChange + ifelse(log2FoldChange > 0, 0.2, -0.2),
              label_y = -log10(pvalue) + 0.5
            )

          p <- p +
            geom_text(
              data = top_genes,
              aes(x = label_x, y = label_y, label = SYMBOL),
              size = label_size, color = txt_col, fontface = "bold",
              check_overlap = TRUE, vjust = 0.5,
              hjust = ifelse(top_genes$log2FoldChange > 0, 0, 1)
            )
        }
      }

      # 保存为SVG
      svg(file, width = 10, height = 8)
      print(p)
      dev.off()
    },
    contentType = "image/svg+xml"
  )

  # 🆕 导出TF靶基因数据（散点图数据）
  output$download_tf_scatter_data <- downloadHandler(
    filename = function() {
      req(selected_tf_targets())
      tf_name <- selected_tf_targets()$tf_name
      paste0("TF_Target_Genes_", tf_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(selected_tf_targets())

      # 获取数据
      data_list <- selected_tf_targets()
      df <- data_list$data
      tf_name <- data_list$tf_name

      # 选择并重命名列
      export_df <- df %>%
        select(
          GeneSymbol = SYMBOL,
          Log2FoldChange = log2FoldChange,
          PValue = pvalue,
          Padj = padj,
          Regulation_Mode = MOR,
          Match_Status = Match_Status,
          Source = source
        ) %>%
        arrange(PValue)

      # 添加元数据
      metadata <- data.frame(
        Parameter = c(
          "TF_Name",
          "N_Targets",
          "N_Consistent",
          "N_Inconsistent",
          "Export_Date"
        ),
        Value = c(
          tf_name,
          nrow(export_df),
          sum(export_df$Match_Status == "Consistent", na.rm = TRUE),
          sum(export_df$Match_Status == "Inconsistent", na.rm = TRUE),
          Sys.Date()
        )
      )

      # 写入元数据和主数据
      write.csv(metadata, file, row.names = FALSE)
      write.table("\n\n=== Target Gene Expression Data ===\n", file, append = TRUE, row.names = FALSE, col.names = FALSE)
      write.csv(export_df, file, append = TRUE, row.names = FALSE)
    },
    contentType = "text/csv"
  )

  # 返回 TF 活性结果供其他模块使用
  return(tf_activity_results)
}