# Metagenomics & Microbiome Pipeline

## 分析流程

### 1. 分类学分析
- **bio-metagenomics-kraken-classification**: Kraken2 快速分类
- **bio-metagenomics-metaphlan-profiling**: MetaPhlAn 物种丰度

### 2. 扩增子分析（16S/ITS）
- **bio-microbiome-amplicon-processing**: DADA2/Deblur 去噪
- **bio-microbiome-taxonomy-assignment**: 物种注释
- **bio-microbiome-qiime2-workflow**: QIIME2 完整流程

### 3. 多样性分析
- **bio-microbiome-diversity-analysis**: Alpha/Beta 多样性

### 4. 差异丰度
- **bio-microbiome-differential-abundance**: DESeq2、ANCOM 差异分析

### 5. 功能预测
- **bio-metagenomics-functional-profiling**: HUMAnN3 功能通路
- **bio-microbiome-functional-prediction**: PICRUSt2 功能预测

### 6. 宏基因组组装
- **bio-genome-assembly-metagenome-assembly**: MetaSPAdes 组装
- **bio-metagenomics-abundance-estimation**: Bin 丰度估计
