  # YuanSeq

<div align="center">

**An Integrated R/Shiny Platform for RNA-seq, Microarray & Single-Cell Downstream Analysis**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E4.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Shiny-green.svg)](https://shiny.posit.co/)

**👨‍💻 Developer**

**Yu Qiao (乔宇)**

PhD in Pharmacology · School of Pharmacy, Shanghai Jiao Tong University

**Supervisors:** [Prof. Feng Qian](https://pharm.sjtu.edu.cn/szdy/2862.html), [Prof. Lei Sun](https://pharm.sjtu.edu.cn/szdy/2870.html)

*YuanSeq is developed as part of the **Yuanclaw** project at SJTU School of Pharmacy*

</div>

  ---

## 📖 Overview

**YuanSeq** is a modular, all-in-one bioinformatics analysis platform designed for comprehensive downstream analysis of bulk RNA-seq, microarray, and **single-cell cluster-derived pseudobulk data**.

### Key Features

| Module | Description | Methods |
|--------|-------------|---------|
| 🔬 **Differential Expression** | RNA-seq / microarray DE analysis | limma-voom, edgeR; 1v1 / nvn comparisons |
| 🧬 **Functional Enrichment** | KEGG, GO, GSEA analysis | clusterProfiler, local KEGG (no rate limit) |
| 🛤️ **Pathway Activity** | Pathway activity inference | ULM, WMEAN, AUCell, GSVA (decoupleR) |
| 🔬 **TF Activity** | Transcription factor activity | CollecTRI network + decoupleR |
| 🆕 **Single-Cell Support** | Downstream analysis for scRNA-seq cluster markers | Use cluster DEGs for enrichment & activity inference |
| 🤖 **AI Interpretation** | AI-powered biological interpretation | Multi-API support (see below) |
| 📊 **Visualization** | Interactive plots with sci-fi themed UI | ggplot2, plotly, pheatmap |

### 🆕 Single-Cell Downstream Analysis

YuanSeq supports downstream analysis of **single-cell RNA-seq data**:
- Import cluster marker genes from Seurat, Scanpy, or other scRNA-seq tools
- Perform KEGG/GO/GSEA enrichment on cluster-specific DEGs
- Infer pathway and TF activity for each cell cluster
- Compare biological programs across cell subpopulations

  ---

## 🤖 AI Interpretation

YuanSeq integrates AI-powered biological interpretation, supporting multiple LLM providers:

| Provider | Models | Features |
|----------|--------|----------|
| **[DeepSeek](https://www.deepseek.com/)** | deepseek-chat, deepseek-reasoner | Cost-effective, excellent Chinese support |
| **[OpenAI](https://openai.com/)** | gpt-4o, gpt-4-turbo, gpt-3.5-turbo | Most capable general models |
| **[Zhipu AI (智谱)](https://www.bigmodel.cn/)** | glm-4, glm-4-flash, glm-4-plus | Chinese NLP leader, fast response |
| **[Local Models](https://ollama.com/)** | Any OpenAI-compatible local model | Data privacy, no external API calls |
| **Custom API** | Any OpenAI-compatible endpoint | Enterprise deployments, self-hosted |

**AI Interpretation Features:**
- Context-aware analysis based on sample metadata (organism, tissue, treatment)
- Multi-format export: Markdown, HTML (with embedded figures), PDF
- Real-time progress indicator
- Conversation history tracking

> ⚠️ **Data Security**: External API calls send data to third-party servers. For sensitive data, use local models or enterprise APIs.

  ### 界面与扩展
  - 科幻主题、玻璃拟态、响应式布局
  - GSEA 模块内提示可配合 [GPSAdb](https://www.gpsadb.com/) fastGPSA 做延伸分析

## 🚀 Installation & Launch

### Requirements

- R >= 4.0
- RStudio (recommended)

### 方式一：GitHub 一键安装 | Method 1: GitHub Install

```r
# 安装包 | Install package
remotes::install_github("Passpoor/Yuanseq")

# 启动应用（自动打开浏览器）| Launch app (auto-opens browser)
YuanSeq::run_app()
```

### 方式二：本地安装 | Method 2: Local Install

```bash
git clone https://github.com/Passpoor/Yuanseq.git
cd Yuanseq
```

```r
# 一键安装依赖 | One-line install dependencies
source("install.R")

# 启动应用 | Launch app
source("run.R")
```

### 手动安装 | Manual Installation

```r
# CRAN packages
install.packages(c(
  "shiny", "shinyjs", "bslib", "ggplot2", "dplyr", "DT",
  "pheatmap", "plotly", "colourpicker", "shinyWidgets", "rlang",
  "tibble", "tidyr", "ggrepel", "RColorBrewer", "VennDiagram",
  "grid", "gridExtra", "httr", "jsonlite", "base64enc"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "AnnotationDbi", "clusterProfiler",
  "org.Mm.eg.db", "org.Hs.eg.db", "GseaVis", "enrichplot",
  "decoupleR", "sva"
))

# Optional: Local KEGG enrichment
remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")

# Launch
shiny::runApp("app.R", launch.browser = TRUE)
```

# Optional: Local KEGG enrichment (recommended)
remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")
```

### Configure AI API (Optional)

```bash
# Windows
mkdir %USERPROFILE%\.yuanseq
copy api_config.example.json %USERPROFILE%\.yuanseq\api_config.json

# Mac/Linux
mkdir -p ~/.yuanseq
cp api_config.example.json ~/.yuanseq/api_config.json
```

Edit `api_config.json`:
```json
{
  "provider": "deepseek",
  "api_key": "YOUR_API_KEY",
  "model": "deepseek-chat"
}
```

### Launch

```r
shiny::runApp("app.R")
```

  # KEGG 本地富集（可选，推荐从 GitHub 安装）
  remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")

## 🙏 Acknowledgments

YuanSeq integrates and acknowledges the following open-source projects:

| Category | Packages | Use |
|----------|----------|-----|
| **Framework** | [shiny](https://cran.r-project.org/package=shiny), [shinyjs](https://cran.r-project.org/package=shinyjs), [bslib](https://cran.r-project.org/package=bslib), [DT](https://cran.r-project.org/package=DT), [plotly](https://cran.r-project.org/package=plotly) | App framework & UI |
| **DE Analysis** | [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/) | Differential expression |
| **Enrichment** | [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [enrichplot](https://bioconductor.org/packages/enrichplot/), [GseaVis](https://bioconductor.org/packages/GseaVis/) | GO/KEGG/GSEA |
| **Local KEGG** | [biofree.qyKEGGtools](https://github.com/Passpoor/biofree.qyKEGGtools) | Offline KEGG enrichment |
| **Activity** | [decoupleR](https://bioconductor.org/packages/decoupleR/) | Pathway/TF activity |
| **Visualization** | [ggplot2](https://cran.r-project.org/package=ggplot2), [pheatmap](https://cran.r-project.org/package=pheatmap), [ggrepel](https://cran.r-project.org/package=ggrepel) | Plotting |
| **AI** | [httr](https://cran.r-project.org/package=httr), [jsonlite](https://cran.r-project.org/package=jsonlite) | API calls |

Thanks to R, Bioconductor, and all package developers!

  ┌──────────┬─────────────┬───────────────┐
  │  提供商  │ provider 值 │   推荐模型    │
  ├──────────┼─────────────┼───────────────┤
  │ DeepSeek │ deepseek    │ deepseek-chat │
  ├──────────┼─────────────┼───────────────┤
  │ OpenAI   │ openai      │ gpt-4o        │
  ├──────────┼─────────────┼───────────────┤
  │ 智谱AI   │ zhipu       │ glm-4-flash   │
  ├──────────┼─────────────┼───────────────┤
  │ 本地模型 │ local       │ custom        │
  └──────────┴─────────────┴───────────────┘

## 📁 Project Structure

```
Yuanseq/
├── app.R                      # Main entry
├── api_config.example.json    # API config template
├── modules/
│   ├── ui_theme.R             # Theme & layout
│   ├── data_input.R           # Data upload
│   ├── differential_analysis.R
│   ├── kegg_enrichment.R
│   ├── gsea_analysis.R
│   ├── pathway_activity.R
│   ├── tf_activity.R
│   ├── ai_interpretation.R    # AI interpretation
│   └── venn_diagram.R
├── workflow/                  # CLI scripts
└── docs/                      # Documentation
```

  shiny::runApp("app.R")
  或使用项目内脚本：launch_app.R、run_app.bat / run_app.sh。

## 📄 License

[MIT License](LICENSE)
