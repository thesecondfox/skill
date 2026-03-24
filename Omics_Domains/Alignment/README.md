# Alignment & Read Processing Pipeline

## 分析流程

### 1. 序列比对
- **bio-read-alignment-bwa-alignment**: BWA 比对（DNA-seq）
- **bio-read-alignment-bowtie2-alignment**: Bowtie2 比对
- **bio-read-alignment-star-alignment**: STAR 比对（RNA-seq）
- **bio-read-alignment-hisat2-alignment**: HISAT2 比对

### 2. 比对文件处理
- **bio-alignment-files-alignment-sorting**: SAM/BAM 排序
- **bio-alignment-files-alignment-indexing**: BAM 索引
- **bio-alignment-files-duplicate-handling**: 去除 PCR 重复
- **bio-alignment-files-bam-statistics**: 比对统计

### 3. 质控
- **bio-read-qc-quality-reports**: FastQC 质量报告
- **bio-read-qc-adapter-trimming**: 接头去除
- **bio-read-qc-quality-filtering**: 质量过滤

### 4. 工具库
- **pysam**: Python SAM/BAM 操作库
