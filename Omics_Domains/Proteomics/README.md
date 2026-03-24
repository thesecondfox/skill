# Proteomics Pipeline

## 分析流程

### 1. 数据导入
- **bio-proteomics-data-import**: 读取 mzML、RAW 等质谱数据

### 2. 肽段鉴定与蛋白推断
- **bio-proteomics-peptide-identification**: 肽段鉴定（搜库）
- **bio-proteomics-protein-inference**: 蛋白推断

### 3. 定量分析
- **bio-proteomics-quantification**: Label-free、TMT、SILAC 定量
- **bio-proteomics-dia-analysis**: DIA 数据分析

### 4. 质控与差异分析
- **bio-proteomics-proteomics-qc**: 质控指标
- **bio-proteomics-differential-abundance**: 差异丰度分析

### 5. 翻译后修饰
- **bio-proteomics-ptm-analysis**: 磷酸化、乙酰化等 PTM 分析

### 6. 工具库
- **pyopenms**: OpenMS Python 接口
