# Spatial Transcriptomics & Spatial Omics Pipeline

## 分析流程

### 1. 数据导入与预处理
- **bio-spatial-transcriptomics-spatial-data-io**: 读取 Visium、MERFISH、seqFISH 等空间组学数据
- **bio-spatial-transcriptomics-spatial-preprocessing**: 质控、归一化、空间坐标处理
- **bio-spatial-transcriptomics-image-analysis**: 组织切片图像分析

### 2. 空间邻域分析
- **bio-spatial-transcriptomics-spatial-neighbors**: 构建空间邻域图
- **bio-spatial-transcriptomics-spatial-statistics**: 空间自相关、Moran's I 等统计

### 3. 空间域识别
- **bio-spatial-transcriptomics-spatial-domains**: 识别空间功能域（SpaGCN、BayesSpace）

### 4. 空间细胞通讯
- **bio-spatial-transcriptomics-spatial-communication**: 空间细胞间通讯分析

### 5. 空间反卷积
- **bio-spatial-transcriptomics-spatial-deconvolution**: Spot 级别的细胞类型反卷积（Cell2location、RCTD）

### 6. 多组学整合
- **bio-spatial-transcriptomics-spatial-multiomics**: 空间多组学数据整合
- **bio-spatial-transcriptomics-spatial-proteomics**: 空间蛋白质组学（CODEX、IMC）

### 7. 可视化
- **bio-spatial-transcriptomics-spatial-visualization**: 空间表达图、热图、交互式可视化

### 8. 专用工具
- **spatial-microbiome**: 空间微生物组分析
- **bio-imaging-mass-cytometry-spatial-analysis**: 成像质谱流式空间分析
