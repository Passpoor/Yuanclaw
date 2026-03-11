# =====================================================
# YuanSeq 一键安装脚本 | One-click Installation Script
# =====================================================
# 运行方式 | Run: source("install.R")
# =====================================================

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║        YuanSeq 安装脚本 | Installation Script         ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# 检查 R 版本 | Check R version
if (R.version.string < "4.0.0") {
  stop("❌ R 版本过低，请升级到 R >= 4.0 | R version too old, please upgrade to R >= 4.0")
}
cat("✅ R 版本检查通过 | R version check passed\n\n")

# 1. 安装 CRAN 包 | Install CRAN packages
cat("📦 安装 CRAN 包 | Installing CRAN packages...\n")
cran_packages <- c(
  "shiny", "shinyjs", "bslib", "ggplot2", "dplyr", "DT",
  "pheatmap", "plotly", "colourpicker", "shinyWidgets", "rlang",
  "tibble", "tidyr", "ggrepel", "RColorBrewer", "VennDiagram",
  "grid", "gridExtra", "httr", "jsonlite", "base64enc", "remotes",
  "later", "sourcetools", "htmltools"
)

missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("   安装 | Installing:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
cat("✅ CRAN 包安装完成 | CRAN packages installed\n\n")

# 2. 安装 BiocManager | Install BiocManager
cat("📦 安装 Bioconductor 包 | Installing Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "edgeR", "limma", "AnnotationDbi", "clusterProfiler",
  "org.Mm.eg.db", "org.Hs.eg.db", "GseaVis", "enrichplot",
  "decoupleR", "sva", "DOSE", "msigdb", "AUCell", "GSVA"
)

missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("   安装 | Installing:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
cat("✅ Bioconductor 包安装完成 | Bioconductor packages installed\n\n")

# 3. 安装 GitHub 包（可选）| Install GitHub packages (optional)
cat("📦 安装 GitHub 包（可选）| Installing GitHub packages (optional)...\n")
tryCatch({
  if (!requireNamespace("biofree.qyKEGGtools", quietly = TRUE)) {
    remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")
    cat("✅ biofree.qyKEGGtools 安装完成 | biofree.qyKEGGtools installed\n")
  } else {
    cat("✅ biofree.qyKEGGtools 已存在 | biofree.qyKEGGtools already installed\n")
  }
}, error = function(e) {
  cat("⚠️ biofree.qyKEGGtools 安装失败，KEGG 本地富集将不可用\n")
  cat("   biofree.qyKEGGtools installation failed, local KEGG enrichment unavailable\n")
})

cat("\n")
cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║            ✅ 安装完成 | Installation Complete          ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("🚀 启动应用 | Launch app:\n")
cat("   source('run.R')\n\n")
