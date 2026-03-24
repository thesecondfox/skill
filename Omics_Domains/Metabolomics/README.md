# Metabolomics Pipeline

## 分析流程

### 1. 数据预处理
- **bio-metabolomics-xcms-preprocessing**: XCMS peak picking
- **bio-metabolomics-msdial-preprocessing**: MS-DIAL 预处理

### 2. 质控与归一化
- **bio-metabolomics-normalization-qc**: 归一化与质控

### 3. 代谢物注释
- **bio-metabolomics-metabolite-annotation**: 数据库匹配注释

### 4. 统计分析
- **bio-metabolomics-statistical-analysis**: 差异代谢物分析

### 5. 通路映射
- **bio-metabolomics-pathway-mapping**: KEGG 通路富集

### 6. 脂质组学
- **bio-metabolomics-lipidomics**: 脂质组学专用分析
