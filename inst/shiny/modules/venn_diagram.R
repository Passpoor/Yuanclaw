# =====================================================
# 韦恩图模块
# =====================================================

venn_diagram_server <- function(input, output, session) {

  # 动态生成集合输入框
  output$venn_inputs <- renderUI({
    n_sets <- input$venn_sets

    tagList(
      lapply(1:n_sets, function(i) {
        fluidRow(
          column(6,
                 textInput(paste0("venn_name_", i),
                           paste0("集合", i, "名称"),
                           value = paste0("Set", i))
          ),
          column(6,
                 textAreaInput(paste0("venn_data_", i),
                               paste0("集合", i, "数据"),
                               placeholder = "每行一个基因/通路，或逗号分隔",
                               rows = 5)
          )
        )
      })
    )
  })

  # 处理韦恩图数据
  venn_data <- eventReactive(input$generate_venn, {
    req(input$venn_sets)

    n_sets <- input$venn_sets
    set_list <- list()
    set_names <- c()

    for (i in 1:n_sets) {
      name_input <- input[[paste0("venn_name_", i)]]
      data_input <- input[[paste0("venn_data_", i)]]

      if (!is.null(name_input) && name_input != "" &&
          !is.null(data_input) && data_input != "") {

        # 处理输入数据：支持换行分隔或逗号分隔
        elements <- unlist(strsplit(data_input, "[\n,]"))
        elements <- trimws(elements)
        elements <- elements[elements != ""]

        if (length(elements) > 0) {
          set_list[[name_input]] <- unique(elements)
          set_names <- c(set_names, name_input)
        }
      }
    }

    if (length(set_list) < 2) {
      showNotification("至少需要2个有效集合才能生成韦恩图", type = "error")
      return(NULL)
    }

    return(list(sets = set_list, names = set_names))
  })

  # 生成韦恩图
  output$venn_plot <- renderPlot({
    req(venn_data())

    data <- venn_data()
    set_list <- data$sets
    set_names <- data$names

    # 颜色设置
    colors <- c(input$venn_color1, input$venn_color2, input$venn_color3,
                input$venn_color4, input$venn_color5)[1:length(set_names)]

    # 根据集合数量选择不同的绘图方式
    if (length(set_names) == 2) {
      # 2集合韦恩图
      grid.newpage()
      venn_plot <- draw.pairwise.venn(
        area1 = length(set_list[[1]]),
        area2 = length(set_list[[2]]),
        cross.area = length(intersect(set_list[[1]], set_list[[2]])),
        category = set_names,
        fill = colors,
        alpha = input$venn_alpha,
        cat.col = colors,
        cat.cex = 1.2,
        cex = 1.5,
        cat.pos = c(0, 0),
        cat.dist = 0.05
      )
    } else if (length(set_names) == 3) {
      # 3集合韦恩图
      grid.newpage()
      venn_plot <- draw.triple.venn(
        area1 = length(set_list[[1]]),
        area2 = length(set_list[[2]]),
        area3 = length(set_list[[3]]),
        n12 = length(intersect(set_list[[1]], set_list[[2]])),
        n23 = length(intersect(set_list[[2]], set_list[[3]])),
        n13 = length(intersect(set_list[[1]], set_list[[3]])),
        n123 = length(Reduce(intersect, set_list)),
        category = set_names,
        fill = colors,
        alpha = input$venn_alpha,
        cat.col = colors,
        cat.cex = 1.2,
        cex = 1.5
      )
    } else if (length(set_names) == 4) {
      # 4集合韦恩图
      grid.newpage()
      venn_plot <- draw.quad.venn(
        area1 = length(set_list[[1]]),
        area2 = length(set_list[[2]]),
        area3 = length(set_list[[3]]),
        area4 = length(set_list[[4]]),
        n12 = length(intersect(set_list[[1]], set_list[[2]])),
        n13 = length(intersect(set_list[[1]], set_list[[3]])),
        n14 = length(intersect(set_list[[1]], set_list[[4]])),
        n23 = length(intersect(set_list[[2]], set_list[[3]])),
        n24 = length(intersect(set_list[[2]], set_list[[4]])),
        n34 = length(intersect(set_list[[3]], set_list[[4]])),
        n123 = length(Reduce(intersect, set_list[1:3])),
        n124 = length(Reduce(intersect, set_list[c(1,2,4)])),
        n134 = length(Reduce(intersect, set_list[c(1,3,4)])),
        n234 = length(Reduce(intersect, set_list[2:4])),
        n1234 = length(Reduce(intersect, set_list)),
        category = set_names,
        fill = colors,
        alpha = input$venn_alpha,
        cat.col = colors,
        cat.cex = 1.1,
        cex = 1.3
      )
    } else if (length(set_names) == 5) {
      # 5集合韦恩图
      grid.newpage()
      venn_plot <- draw.quintuple.venn(
        area1 = length(set_list[[1]]),
        area2 = length(set_list[[2]]),
        area3 = length(set_list[[3]]),
        area4 = length(set_list[[4]]),
        area5 = length(set_list[[5]]),
        n12 = length(intersect(set_list[[1]], set_list[[2]])),
        n13 = length(intersect(set_list[[1]], set_list[[3]])),
        n14 = length(intersect(set_list[[1]], set_list[[4]])),
        n15 = length(intersect(set_list[[1]], set_list[[5]])),
        n23 = length(intersect(set_list[[2]], set_list[[3]])),
        n24 = length(intersect(set_list[[2]], set_list[[4]])),
        n25 = length(intersect(set_list[[2]], set_list[[5]])),
        n34 = length(intersect(set_list[[3]], set_list[[4]])),
        n35 = length(intersect(set_list[[3]], set_list[[5]])),
        n45 = length(intersect(set_list[[4]], set_list[[5]])),
        n123 = length(Reduce(intersect, set_list[1:3])),
        n124 = length(Reduce(intersect, set_list[c(1,2,4)])),
        n125 = length(Reduce(intersect, set_list[c(1,2,5)])),
        n134 = length(Reduce(intersect, set_list[c(1,3,4)])),
        n135 = length(Reduce(intersect, set_list[c(1,3,5)])),
        n145 = length(Reduce(intersect, set_list[c(1,4,5)])),
        n234 = length(Reduce(intersect, set_list[2:4])),
        n235 = length(Reduce(intersect, set_list[c(2,3,5)])),
        n245 = length(Reduce(intersect, set_list[c(2,4,5)])),
        n345 = length(Reduce(intersect, set_list[3:5])),
        n1234 = length(Reduce(intersect, set_list[1:4])),
        n1235 = length(Reduce(intersect, set_list[c(1,2,3,5)])),
        n1245 = length(Reduce(intersect, set_list[c(1,2,4,5)])),
        n1345 = length(Reduce(intersect, set_list[c(1,3,4,5)])),
        n2345 = length(Reduce(intersect, set_list[2:5])),
        n12345 = length(Reduce(intersect, set_list)),
        category = set_names,
        fill = colors,
        alpha = input$venn_alpha,
        cat.col = colors,
        cat.cex = 1.0,
        cex = 1.2
      )
    }

    # 存储交集信息用于交互
    venn_intersections <<- calculate_intersections(set_list)

    return(venn_plot)
  })

  # 计算所有可能的交集
  calculate_intersections <- function(set_list) {
    set_names <- names(set_list)
    intersections <- list()

    # 生成所有可能的组合
    for (i in 1:length(set_names)) {
      combos <- combn(set_names, i, simplify = FALSE)
      for (combo in combos) {
        combo_name <- paste(combo, collapse = " & ")
        if (length(combo) == 1) {
          # 单个集合
          elements <- set_list[[combo]]
        } else {
          # 多个集合的交集
          elements <- Reduce(intersect, set_list[combo])
        }
        intersections[[combo_name]] <- elements
      }
    }

    return(intersections)
  }

  # 处理韦恩图点击事件
  observeEvent(input$venn_click, {
    req(venn_intersections)

    # 这里简化处理，实际应用中可以根据点击坐标确定具体区域
    # 现在显示所有交集供用户选择
    showModal(modalDialog(
      title = "选择交集区域",
      selectInput("selected_intersection", "选择交集:",
                  choices = names(venn_intersections)),
      footer = tagList(
        actionButton("confirm_intersection", "确认选择"),
        modalButton("取消")
      )
    ))
  })

  # 确认选择的交集
  observeEvent(input$confirm_intersection, {
    req(input$selected_intersection, venn_intersections)

    selected_elements <- venn_intersections[[input$selected_intersection]]

    # 复制到剪贴板
    if (length(selected_elements) > 0) {
      clip_text <- paste(selected_elements, collapse = "\n")
      writeClipboard(clip_text)
      showNotification(paste("已复制", length(selected_elements), "个元素到剪贴板"),
                       type = "message")
    } else {
      showNotification("选择的交集为空", type = "warning")
    }

    removeModal()
  })

  # 显示交集结果
  output$venn_result_ui <- renderUI({
    req(venn_data())

    data <- venn_data()
    set_list <- data$sets
    set_names <- data$names

    # 计算基本统计信息
    total_elements <- length(unique(unlist(set_list)))
    intersection_counts <- list()

    # 计算所有交集
    for (i in 1:length(set_names)) {
      combos <- combn(set_names, i, simplify = FALSE)
      for (combo in combos) {
        combo_name <- paste(combo, collapse = " ∩ ")
        if (length(combo) == 1) {
          count <- length(set_list[[combo]])
        } else {
          count <- length(Reduce(intersect, set_list[combo]))
        }
        intersection_counts[[combo_name]] <- count
      }
    }

    # 创建结果展示
    tagList(
      h4("韦恩图统计信息"),
      div(class = "venn-result-box",
          p(strong("总唯一元素数: "), total_elements),
          lapply(names(intersection_counts), function(name) {
            div(
              p(strong(name), ": ", intersection_counts[[name]], "个元素",
                actionButton(paste0("copy_", gsub(" ", "_", name)), "复制",
                             class = "btn-xs copy-btn",
                             onclick = sprintf('Shiny.setInputValue("copy_intersection", "%s", {priority: "event"})', name))
              )
            )
          })
      )
    )
  })

  # 处理复制按钮点击
  observeEvent(input$copy_intersection, {
    req(input$copy_intersection, venn_intersections)

    intersection_name <- input$copy_intersection
    selected_elements <- venn_intersections[[intersection_name]]

    if (length(selected_elements) > 0) {
      clip_text <- paste(selected_elements, collapse = "\n")
      writeClipboard(clip_text)
      showNotification(paste("已复制", length(selected_elements), "个元素到剪贴板"),
                       type = "message")
    } else {
      showNotification("选择的交集为空", type = "warning")
    }
  })

  # 显示详细交集数据表
  output$venn_table <- DT::renderDataTable({
    req(venn_intersections)

    # 准备数据框
    intersection_df <- data.frame(
      交集区域 = names(venn_intersections),
      元素数量 = sapply(venn_intersections, length),
      元素列表 = sapply(venn_intersections, function(x) paste(x, collapse = ", ")),
      stringsAsFactors = FALSE
    )

    DT::datatable(intersection_df,
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE) %>%
      formatStyle(columns = 1:3, fontSize = '14px')
  })
}