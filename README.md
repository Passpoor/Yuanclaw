# YuanSeq

<div align="center">

**基于 R/Shiny 的综合生物信息学分析平台**

**A comprehensive bioinformatics analysis platform built with R/Shiny**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E4.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Shiny-green.svg)](https://shiny.posit.co/)

</div>

---

## 📖 简介 | Introduction

**YuanSeq（源Seq）** 是一个模块化的生物信息学分析平台，提供从差异表达分析、功能富集到通路活性推断的完整流程。

**核心能力：**
- 🔬 **差异表达分析**：limma-voom、edgeR；支持 1v1 / nvn 比较
- 🧬 **富集分析**：KEGG（含本地/背景基因）、GO、GSEA
- 🛤️ **通路活性推断**：ULM/WMEAN/AUCell/GSVA（decoupleR）
- 🔬 **转录因子活性**：CollecTRI 网络与 decoupleR
- 🤖 **AI 智能解读**：多 API 支持，生成专业生物学解读报告
- 📊 **交互式可视化**：科幻主题 UI、玻璃拟态设计、响应式布局

---

## 🚀 快速开始 | Quick Start

### 环境要求

- R >= 4.0
- RStudio（推荐）

### 安装步骤

#### 1️⃣ 克隆仓库

```bash
git clone https://github.com/Passpoor/Yuanseq.git
cd Yuanseq
```

#### 2️⃣ 安装依赖包

在 R 中执行：

```r
# CRAN 包
install.packages(c(
  "shiny", "shinyjs", "bslib", "ggplot2", "dplyr", "DT",
  "pheatmap", "plotly", "colourpicker", "shinyWidgets", "rlang",
  "tibble", "tidyr", "ggrepel", "RColorBrewer", "VennDiagram",
  "grid", "gridExtra", "httr", "jsonlite", "base64enc"
))

# Bioconductor 包
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "AnnotationDbi", "clusterProfiler",
  "org.Mm.eg.db", "org.Hs.eg.db", "GseaVis", "enrichplot",
  "decoupleR", "sva"
))

# KEGG 本地富集（可选）
remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")
```

#### 3️⃣ 配置 AI API（可选）

```bash
# Windows
copy api_config.example.json %USERPROFILE%\.yuanseq\api_config.json

# Mac/Linux
mkdir -p ~/.yuanseq
cp api_config.example.json ~/.yuanseq/api_config.json
```

编辑配置文件，填入您的 API Key：

```json
{
  "provider": "deepseek",
  "api_key": "您的API Key",
  "model": "deepseek-chat"
}
```

**支持的 API 提供商：**

| 提供商 | provider | 推荐模型 | 特点 |
|--------|----------|----------|------|
| DeepSeek | `deepseek` | deepseek-chat | 性价比高，中文友好 |
| OpenAI | `openai` | gpt-4o | 综合能力强 |
| 智谱AI | `zhipu` | glm-4-flash | 国产，响应快 |
| 本地模型 | `local` | custom | 数据不出本地 |

> ⚠️ **数据安全**：使用外部 API 会将数据发送到第三方服务器，敏感数据请使用本地模型。

#### 4️⃣ 启动应用

```r
shiny::runApp("app.R")
```

---

## 📋 功能详情 | Features

### 核心分析模块

| 模块 | 功能 | 方法 |
|------|------|------|
| 差异表达 | RNA-seq / 芯片差异分析 | limma-voom, edgeR |
| KEGG 富集 | 通路富集分析 | clusterProfiler, 本地KEGG |
| GO 富集 | 基因本体富集 | clusterProfiler |
| GSEA | 基因集富集分析 | GseaVis, enrichplot |
| 通路活性 | 通路活性推断 | decoupleR |
| TF 活性 | 转录因子活性 | decoupleR + CollecTRI |
| 韦恩图 | 多组交集可视化 | VennDiagram |

### 🆕 AI 解读功能

- **多 API 支持**：OpenAI、DeepSeek、智谱AI、本地模型
- **智能解读**：基于分析结果生成专业生物学报告
- **样本信息**：输入物种、组织、处理条件等，生成针对性解读
- **多格式导出**：Markdown、HTML（含嵌入图片）、PDF
- **实时进度**：显示分析进度条

---

## 🙏 饮水思源 | Acknowledgments

YuanSeq 是集成平台，依赖并致谢以下开源项目：

| 类别 | 包名 | 用途 |
|------|------|------|
| **框架** | [shiny](https://cran.r-project.org/package=shiny), [shinyjs](https://cran.r-project.org/package=shinyjs), [bslib](https://cran.r-project.org/package=bslib) | 应用框架 |
| **差异分析** | [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/) | 差异表达 |
| **富集分析** | [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [enrichplot](https://bioconductor.org/packages/enrichplot/), [GseaVis](https://bioconductor.org/packages/GseaVis/) | GO/KEGG/GSEA |
| **活性推断** | [decoupleR](https://bioconductor.org/packages/decoupleR/) | 通路/TF 活性 |
| **可视化** | [ggplot2](https://cran.r-project.org/package=ggplot2), [pheatmap](https://cran.r-project.org/package=pheatmap), [plotly](https://cran.r-project.org/package=plotly) | 图表绑制 |
| **AI 功能** | [httr](https://cran.r-project.org/package=httr), [jsonlite](https://cran.r-project.org/package=jsonlite) | API 调用 |

感谢 R、Bioconductor 社区及所有上游开发者！

---

## 📁 项目结构 | Structure

```
Yuanseq/
├── app.R                      # 主入口
├── api_config.example.json    # API 配置模板
├── modules/
│   ├── ui_theme.R             # 主题与布局
│   ├── data_input.R           # 数据上传
│   ├── differential_analysis.R
│   ├── kegg_enrichment.R
│   ├── gsea_analysis.R
│   ├── pathway_activity.R
│   ├── tf_activity.R
│   ├── ai_interpretation.R    # AI 解读
│   └── venn_diagram.R
├── workflow/                  # 命令行脚本
└── docs/                      # 文档
```

---

## 👨‍💻 开发者 | Developer

**乔宇 Yu Qiao**

上海交通大学药学院 · 药理学博士

School of Pharmacy, Shanghai Jiao Tong University · PhD in Pharmacology

**导师 Supervisors：**

- [钱峰教授 Prof. Feng Qian](https://pharm.sjtu.edu.cn/szdy/2862.html)
- [孙磊教授 Prof. Lei Sun](https://pharm.sjtu.edu.cn/szdy/2870.html)

---

## 📄 许可证 | License

[MIT License](LICENSE)
