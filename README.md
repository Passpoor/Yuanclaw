# YuanSeq

YuanSeq is a web-based R/Shiny platform for comprehensive bioinformatics analysis of RNA-seq and microarray data: differential expression (limma-voom / edgeR), functional enrichment (KEGG / GO / GSEA), transcription factor and pathway activity inference, and interactive visualization. Developed at Shanghai Jiao Tong University School of Pharmacy.

**Repository:** [https://github.com/Passpoor/Xseq0.1](https://github.com/Passpoor/Xseq0.1)

**开发者 Developer:** 乔宇 Yu Qiao · 上海交通大学药学院 药理学博士 | School of Pharmacy, Shanghai Jiao Tong University · PhD in Pharmacology  

**导师 Supervisors:** [钱峰教授 Prof. Feng Qian](https://pharm.sjtu.edu.cn/szdy/2862.html)、[孙磊教授 Prof. Lei Sun](https://pharm.sjtu.edu.cn/szdy/2870.html)

---

## 项目概述 | About

YuanSeq（源Seq）为模块化生物信息学分析平台，基于 Shiny 开发，提供从差异表达、富集分析到通路活性推断的完整流程，支持科幻主题 UI 与日夜模式切换。本项目集成 R/Bioconductor 社区开源包，饮水思源，在此致谢所有上游开发者。

---

## 功能特性 | Features

### 核心功能
- **差异表达分析**: limma-voom、edgeR；支持 1v1 / nvn 比较
- **富集分析**: KEGG（含本地/背景基因）、GO、GSEA（含 Leading Edge 与 GPSAdb 延伸提示）
- **通路活性推断**: ULM/WMEAN/AUCell/GSVA（decoupleR），基于 KEGG 富集结果
- **转录因子活性**: CollecTRI 网络与 decoupleR
- **韦恩图、火山图**: 多组交集、多种差异结果格式

### 界面与扩展
- 科幻主题、玻璃拟态、响应式布局
- GSEA 模块内提示可配合 [GPSAdb](https://www.gpsadb.com/) fastGPSA 做延伸分析

---

## 安装与运行 | Install & Run

### 要求
- R >= 4.0
- 需安装 Shiny、BiocManager 及以下依赖

### 1. 克隆仓库
```bash
git clone https://github.com/Passpoor/Xseq0.1.git
cd Xseq0.1
```

### 2. 安装 R 包
在 R 中执行：
```r
install.packages(c("shiny", "shinyjs", "bslib", "ggplot2", "dplyr", "DT",
  "pheatmap", "plotly", "colourpicker", "shinyWidgets", "rlang",
  "tibble", "tidyr", "ggrepel", "RColorBrewer", "VennDiagram", "grid", "gridExtra"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "AnnotationDbi", "clusterProfiler",
  "org.Mm.eg.db", "org.Hs.eg.db", "GseaVis", "enrichplot", "decoupleR", "sva"))

# KEGG 本地富集（可选，推荐从 GitHub 安装）
remotes::install_github("Passpoor/biofree.qyKEGGtools", upgrade = "never")
```

### 3. 启动应用
```r
shiny::runApp("app.R")
```
或使用项目内脚本：`launch_app.R`、`run_app.bat` / `run_app.sh`。

---

## 饮水思源 · 致谢 | Acknowledgments

YuanSeq 为集成平台，未重复造轮子，依赖并致谢以下 R/Bioconductor 开源包及社区。

| 类别 | 包名 | 用途 |
|------|------|------|
| **框架与 UI** | [shiny](https://cran.r-project.org/package=shiny), [shinyjs](https://cran.r-project.org/package=shinyjs), [bslib](https://cran.r-project.org/package=bslib), [DT](https://cran.r-project.org/package=DT), [plotly](https://cran.r-project.org/package=plotly), [colourpicker](https://cran.r-project.org/package=colourpicker), [shinyWidgets](https://cran.r-project.org/package=shinyWidgets) | 应用框架与交互界面 |
| **差异分析** | [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/) | RNA-seq / 芯片差异表达 |
| **注释与富集** | [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/), [org.Mm.eg.db](https://bioconductor.org/packages/org.Mm.eg.db/), [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/), [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [enrichplot](https://bioconductor.org/packages/enrichplot/), [GseaVis](https://bioconductor.org/packages/GseaVis/) | 基因注释、GO/KEGG/GSEA 富集与可视化 |
| **KEGG 本地** | [biofree.qyKEGGtools](https://github.com/Passpoor/biofree.qyKEGGtools) | 本地 KEGG 富集（可选） |
| **通路与 TF** | [decoupleR](https://bioconductor.org/packages/decoupleR/) | 通路活性、转录因子活性推断 |
| **可视化** | [ggplot2](https://cran.r-project.org/package=ggplot2), [pheatmap](https://cran.r-project.org/package=pheatmap), [ggrepel](https://cran.r-project.org/package=ggrepel), [RColorBrewer](https://cran.r-project.org/package=RColorBrewer), [VennDiagram](https://cran.r-project.org/package=VennDiagram), [grid](https://cran.r-project.org/package=grid), [gridExtra](https://cran.r-project.org/package=gridExtra) | 图表与排版 |
| **数据处理** | [dplyr](https://cran.r-project.org/package=dplyr), [tibble](https://cran.r-project.org/package=tibble), [tidyr](https://cran.r-project.org/package=tidyr), [rlang](https://cran.r-project.org/package=rlang), [later](https://cran.r-project.org/package=later) | 数据整理与异步 |

芯片分析模块另用 [reshape2](https://cran.r-project.org/package=reshape2)、[sva](https://bioconductor.org/packages/sva/) 等。

感谢 R、Bioconductor 及上述所有包的开发者与维护者。

---

## 项目结构 | Structure

```
├── app.R                 # 主入口
├── config/               # 配置
├── modules/              # Shiny 模块
│   ├── ui_theme.R        # 主题与布局
│   ├── data_input.R      # 数据上传与注释
│   ├── differential_analysis.R
│   ├── kegg_enrichment.R
│   ├── gsea_analysis.R
│   ├── pathway_activity.R # 通路活性推断
│   ├── tf_activity.R
│   └── venn_diagram.R
├── workflow/             # 工作流脚本
├── tests/                # 测试
└── docs/                 # 文档
```

---

## 许可证 | License

MIT License. See [LICENSE](LICENSE) for details.
