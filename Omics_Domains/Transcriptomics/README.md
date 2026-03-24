# Transcriptomics (RNA-seq) Pipeline

## 分析流程

### 1. 数据质控
- **bio-read-qc-quality-reports**: FastQC 质量报告
- **bio-read-qc-adapter-trimming**: 接头去除（Trimmomatic、Cutadapt）
- **bio-read-qc-rnaseq-qc**: RNA-seq 特异性质控

### 2. 序列比对
- **bio-read-alignment-star-alignment**: STAR 比对（推荐）
- **bio-read-alignment-hisat2-alignment**: HISAT2 比对

### 3. 定量
- **bio-rna-quantification-alignment-free-quant**: Salmon/Kallisto 无比对定量（快速）
- **bio-rna-quantification-featurecounts-counting**: 基于比对的计数
- **bio-rna-quantification-tximport-workflow**: 导入定量结果

### 4. 差异表达分析
- **bio-differential-expression-deseq2-basics**: DESeq2 分析（推荐）
- **bio-differential-expression-edger-basics**: edgeR 分析
- **bio-differential-expression-batch-correction**: 批次效应校正
- **bio-differential-expression-de-visualization**: 火山图、MA plot 可视化
- **pydeseq2**: Python 版 DESeq2

### 5. 可变剪接分析
- **bio-alternative-splicing-splicing-quantification**: 剪接事件定量（rMATS、SUPPA）
- **bio-alternative-splicing-differential-splicing**: 差异剪接分析
- **bio-alternative-splicing-isoform-switching**: 亚型切换分析

### 6. 其他 RNA-seq 类型
- **bio-small-rna-seq-smrna-preprocessing**: 小 RNA 预处理
- **bio-small-rna-seq-differential-mirna**: miRNA 差异分析
- **bio-ribo-seq-riboseq-preprocessing**: Ribo-seq 核糖体足迹分析
- **bio-ribo-seq-translation-efficiency**: 翻译效率计算
