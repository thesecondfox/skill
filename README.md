# Bioinformatics Skills Collection

492 curated skills for computational biology and bioinformatics, organized in a Two-Tier architecture.

## Architecture

```
├── Common_Skills/       (248)  跨领域通用工具
└── Omics_Domains/       (244)  按组学领域分类，每个子目录含 Pipeline README
```

## Omics Domains (17 领域)

| 领域 | Skills | 覆盖范围 |
|------|--------|---------|
| Genomics | 39 | 变异检测、基因组组装与注释、群体遗传学、结构变异 |
| Epigenomics | 26 | ChIP-seq、ATAC-seq、DNA 甲基化、Hi-C、衰老时钟 |
| Transcriptomics | 24 | RNA-seq、差异表达、可变剪接、Ribo-seq、小 RNA |
| Immunology | 17 | TCR/BCR 分析、流式细胞术、免疫信息学、新抗原预测 |
| Single_Cell | 16 | 预处理、聚类、注释、轨迹推断、扰动预测（27 模型基准） |
| Spatial | 15 | 空间转录组、空间蛋白组、空间域检测、反卷积 |
| Alignment | 14 | BWA/STAR/HISAT2 比对、SAM/BAM 操作、多序列比对 |
| Functional_Genomics | 14 | CRISPR 筛选、CLIP-seq、Perturb-seq |
| Metagenomics | 14 | 物种分类、功能注释、耐药基因检测、菌株追踪 |
| Clinical | 10 | 临床变异解读、液体活检、药物基因组学、肿瘤突变负荷 |
| Pathway_Analysis | 10 | GO/KEGG/Reactome 富集、GSEA、WikiPathways |
| Proteomics | 10 | DIA 分析、蛋白定量、翻译后修饰、差异丰度 |
| Chemoinformatics | 9 | 分子描述符、ADMET 预测、虚拟筛选、相似性搜索 |
| Machine_Learning | 8 | 组学分类器、生物标志物发现、生存分析、模型验证 |
| Metabolomics | 8 | XCMS/MS-DIAL 预处理、代谢物注释、脂质组学 |
| Structural_Biology | 6 | AlphaFold 预测、PDB 操作、分子动力学 |
| Long_Read | 4 | Nanopore/HiFi 比对、Clair3 变异检测、甲基化检测 |

每个领域目录下都有 README.md，按分析 Pipeline 的生命周期串联所有 skill。

## Common Skills (248) 分类概览

### 数据库访问
PubMed, UniProt, Ensembl, KEGG, Reactome, STRING, PDB, ClinVar, gnomAD, COSMIC, DrugBank, ChEMBL, PubChem, GEO, GTEx, GWAS Catalog, DepMap, OpenTargets 等 30+ 数据库

### 序列与基因组操作
序列读写与格式转换、基因组区间运算（BED/GTF/BigWig）、表达矩阵处理、比较基因组学、限制性酶切分析

### 可视化
Circos、基因组浏览器、热图、火山图、UpSet 图、网络图、多面板组图、交互式可视化、配色方案

### 实验设计与统计
样本量计算、功效分析、批次设计、多重检验校正、统计分析、探索性数据分析

### 工作流管理
Snakemake、Nextflow、CWL、WDL 及 40+ 预构建 Pipeline（RNA-seq、ATAC-seq、GWAS、宏基因组等）

### 科研写作与报告
科学写作、文献管理、LaTeX/PPTX 海报、Jupyter/Quarto/RMarkdown 报告、图表导出

### 核心库
scanpy, anndata, biopython, pysam, scikit-bio, rdkit, deepchem, pydeseq2, scvi-tools, scvelo, cobrapy, networkx, matplotlib, seaborn, plotly, scikit-learn, statsmodels, polars, dask

### 专题分析
表观转录组学（m6A）、生态基因组学、流行病基因组学、基因组工程（CRISPR 设计）、系统生物学（代谢建模）、时间序列基因组学

## Quick Start

1. 确定分析类型，进入对应的 Omics_Domains 子目录
2. 阅读该领域的 README.md，了解完整 Pipeline
3. 根据当前分析阶段选择对应 skill

## Changelog

详见 [CHANGELOG_20260326.md](CHANGELOG_20260326.md)
