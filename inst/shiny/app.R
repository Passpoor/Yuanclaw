# =====================================================
# YuanSeq
# 开发者 Developer: 乔宇 Yu Qiao
# 上海交通大学药学院 药理学博士
# PhD in Pharmacology, School of Pharmacy, Shanghai Jiao Tong University
# 导师 Supervisors: 钱峰教授 Prof. Feng Qian、孙磊教授 Prof. Lei Sun
# =====================================================

# 设置上传大小限制
options(shiny.maxRequestSize = 100 * 1024^2) # 将上传上限设置为 100MB

# 加载必要的包
library(shiny)
library(shinyjs)
library(bslib)
library(ggplot2)
library(dplyr)
library(DT)
library(pheatmap)
library(plotly)
library(colourpicker)
library(shinyWidgets)
library(rlang)
library(later)

# 生物信包加载
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(clusterProfiler)
  try(library(org.Mm.eg.db), silent=TRUE)
  try(library(org.Hs.eg.db), silent=TRUE)
  try(library(biofree.qyKEGGtools), silent=TRUE)
  try(library(GseaVis), silent=TRUE)
  try(library(enrichplot), silent=TRUE)

  # === decoupleR 模块所需包 ===
  library(decoupleR)
  library(tibble)
  library(tidyr)
  library(ggrepel)
  library(RColorBrewer)

  # === 韦恩图所需包 ===
  library(VennDiagram)
  library(grid)
  library(gridExtra)
})

# ===============================
# 加载模块
# =====================================================

# 加载配置
source("config/config.R")

# #region agent log (debug)
debug_log_ndjson <- function(location, message, data = list()) {
  try({
    line <- paste0(
      "{\"sessionId\":\"aec83c\",\"timestamp\":", as.integer(Sys.time()) * 1000,
      ",\"location\":\"", location,
      "\",\"message\":\"", gsub("\"", "'", message),
      "\",\"data\":", paste0("{", paste0(sprintf("\"%s\":\"%s\"", names(data), as.character(data)), collapse = ","), "}"),
      "}\n"
    )
    cat(line, file = "debug-aec83c.log", append = TRUE)
  }, silent = TRUE)
}
# #endregion

# 加载核心模块
tryCatch({
  source("modules/ui_theme.R")
  debug_log_ndjson("app.R:source(ui_theme)", "ui_theme_loaded", list(ok = "true"))
}, error = function(e) {
  debug_log_ndjson("app.R:source(ui_theme)", "ui_theme_source_error", list(error = conditionMessage(e)))
  stop(e)
})
source("modules/data_input.R")
source("modules/differential_analysis.R")
source("modules/kegg_enrichment.R")
source("modules/go_analysis.R")   # GO分析模块
source("modules/gsea_analysis.R")
source("modules/tf_activity.R")
source("modules/pathway_activity.R")  # 🆕 通路活性分析模块
source("modules/chip_analysis.R")      # 🆕 芯片数据分析模块
source("modules/venn_diagram.R")
source("modules/ai_interpretation.R")   # 🤖 AI 解读模块

# ===============================
# 主应用
# =====================================================

# 创建UI
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    sci_fi_css,
    tags$style(HTML("
      body { color: inherit; }
      .small-box { color: #fff !important; }
      .shiny-notification {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        border-radius: 10px;
        backdrop-filter: blur(10px);
      }
    "))
  ),
  uiOutput("app_ui")
)

# 创建Server
server <- function(input, output, session) {

  # 设置初始主题 - 使用默认主题，依赖CSS夜间模式
  initial_theme <- bs_theme(version = 5)

  # 动态渲染 UI
  output$app_ui <- renderUI({
    main_app_ui(initial_theme)
  })

  # =====================================================
  # 主题切换逻辑
  # =====================================================

  observeEvent(input$theme_toggle, {
    if(input$theme_toggle) {
      # 夜间模式
      session$sendCustomMessage("toggle-darkmode", TRUE)
    } else {
      # 日间模式
      session$sendCustomMessage("toggle-darkmode", FALSE)
    }
  }, ignoreInit = TRUE)

  # =====================================================
  # 调用各功能模块
  # =====================================================

  # 数据输入模块
  data_input_server(input, output, session)

  # 差异分析模块
  deg_results <- differential_analysis_server(input, output, session)

  # KEGG富集模块
  kegg_results <- kegg_enrichment_server(input, output, session, deg_results)

  # GO富集分析模块
  go_results <- go_analysis_server(input, output, session, deg_results)

  # GSEA分析模块
  gsea_analysis_server(input, output, session, deg_results)

  # 转录因子活性模块
  tf_activity_results <- tf_activity_server(input, output, session, deg_results)

  # 🆕 通路活性分析模块
  pathway_activity_results <- pathway_activity_server(input, output, session, deg_results, kegg_results)

  # 韦恩图模块
  venn_diagram_server(input, output, session)

  # 🆕 芯片数据分析模块
  chip_analysis_server(input, output, session, deg_results)

  # 🤖 AI 解读模块
  ai_interpretation_server(
    input, output, session,
    deg_results = deg_results,
    kegg_results = kegg_results,
    go_results = go_results,
    tf_activity_results = tf_activity_results,
    pathway_activity_results = pathway_activity_results
  )

}

# =====================================================
# 🚀 启动应用
# =====================================================
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))