# =====================================================
# YuanSeq 一键启动脚本 | One-click Launch Script
# =====================================================
# 运行方式 | Run: source("run.R")
# =====================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║              YuanSeq 启动中...                        ║\n")
cat("║              Launching YuanSeq...                     ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# 获取应用路径 | Get app path
app_path <- getwd()

# 启动 Shiny 应用，自动打开浏览器 | Launch Shiny app with browser
shiny::runApp(
  appDir = app_path,
  launch.browser = TRUE,  # 自动打开浏览器 | Auto-open browser
  port = 3838,            # 固定端口 | Fixed port
  host = "127.0.0.1"
)
