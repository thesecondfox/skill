# 生物信息学专用 Skills 目录

## 📊 筛选统计

- **原始总数**: 604 个
- **第一轮去重**: 569 个（精简 35 个）
- **第二轮生信筛选**: 490 个（移除 79 个非生信）
- **最终保留**: 490 个纯生信 skills

## 📁 目录结构

```
skills_bioinformatics_only/
├── Common_Skills/          248 个生信通用工具
│   ├── 生信数据库（PubMed, UniProt, KEGG, GnomAD 等）
│   ├── 生信可视化（matplotlib, seaborn, plotly）
│   ├── 序列操作（biopython, scikit-bio）
│   ├── 统计分析（scikit-learn, statsmodels）
│   └── 科研写作（LaTeX, 文献管理）
│
└── Omics_Domains/          242 个组学专业工具
    ├── Single_Cell/        15 个
    ├── Spatial/            15 个
    ├── Transcriptomics/    24 个
    ├── Genomics/           39 个
    ├── Epigenomics/        25 个
    ├── Proteomics/         10 个
    ├── Metabolomics/       8 个
    ├── Metagenomics/       14 个
    └── ... 其他 9 个领域
```

## 🎯 移除的非生信内容

已移除以下类别（79 个）：
- 金融数据工具（alpha-vantage, edgartools, fred-economic-data）
- 通用 AI 工具（consciousness-council, hypothesis-generation）
- 量子计算（qiskit, cirq, pennylane, qutip）
- 通用编程（modal, simpy, pufferlib）
- 非生信数据源（datacommons, usfiscaldata, hedgefundmonitor）

## ✅ 保留的生信核心内容

### 组学分析全覆盖
- 单细胞、空间、转录组、基因组、表观、蛋白、代谢、宏基因组
- 免疫组库、功能基因组、结构生物学、化学信息学

### 生信数据库完整
- 序列数据库：PubMed, UniProt, Ensembl, PDB
- 变异数据库：ClinVar, gnomAD, COSMIC
- 通路数据库：KEGG, Reactome, STRING
- 药物数据库：DrugBank, ChEMBL, PubChem

### 科研工作流
- 科学写作、可视化、文献管理
- 统计分析、机器学习
- 工作流管理（Snakemake, Nextflow）

## 📖 使用方式

新目录位置：`/home/yzhou/.claude/skills_bioinformatics_only/`

每个组学子目录都有 README.md，按 Pipeline 流程组织。
