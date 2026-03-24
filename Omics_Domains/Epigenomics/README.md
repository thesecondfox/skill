# Epigenomics Pipeline

## 分析流程

### 1. ChIP-seq 分析
- **bio-chip-seq-chipseq-qc**: ChIP-seq 质控
- **bio-chip-seq-peak-calling**: Peak calling（MACS2）
- **bio-chip-seq-peak-annotation**: Peak 注释到基因
- **bio-chip-seq-differential-binding**: 差异结合分析
- **bio-chip-seq-motif-analysis**: Motif 富集分析

### 2. ATAC-seq 分析
- **bio-atac-seq-atac-qc**: ATAC-seq 质控
- **bio-atac-seq-atac-peak-calling**: 开放染色质 peak calling
- **bio-atac-seq-differential-accessibility**: 差异可及性分析
- **bio-atac-seq-footprinting**: 转录因子足迹分析

### 3. DNA 甲基化分析
- **bio-methylation-analysis-bismark-alignment**: Bisulfite-seq 比对
- **bio-methylation-analysis-methylation-calling**: 甲基化位点检测
- **bio-methylation-analysis-dmr-detection**: 差异甲基化区域（DMR）
- **bio-methylation-analysis-methylkit-analysis**: MethylKit 分析流程

### 4. Hi-C 染色质互作
- **bio-hi-c-analysis-hic-data-io**: Hi-C 数据读取
- **bio-hi-c-analysis-tad-detection**: TAD 检测
- **bio-hi-c-analysis-loop-calling**: 染色质环检测
- **bio-hi-c-analysis-compartment-analysis**: A/B compartment 分析

### 5. 工具库
- **deeptools**: ChIP-seq/ATAC-seq 可视化工具集
