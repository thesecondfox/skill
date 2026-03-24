# Functional Genomics Pipeline

## 分析流程

### 1. CRISPR 筛选
- **bio-crispr-screens-screen-qc**: 筛选质控
- **bio-crispr-screens-mageck-analysis**: MAGeCK 分析
- **bio-crispr-screens-hit-calling**: Hit 基因识别

### 2. CLIP-seq（RNA-蛋白互作）
- **bio-clip-seq-clip-preprocessing**: CLIP 预处理
- **bio-clip-seq-clip-peak-calling**: 结合位点检测
- **bio-clip-seq-clip-motif-analysis**: Motif 分析

### 3. 扰动实验
- **bio-single-cell-perturb-seq**: Perturb-seq 分析
