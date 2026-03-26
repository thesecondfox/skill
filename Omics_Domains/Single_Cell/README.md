# Single Cell Analysis Pipeline

## 分析流程

### 1. 数据导入与预处理
- **bio-single-cell-data-io**: 读取 10X、H5AD 等格式的单细胞数据
- **bio-single-cell-preprocessing**: 质控、过滤低质量细胞、归一化
- **bio-single-cell-doublet-detection**: 检测和去除双细胞

### 2. 降维与聚类
- **bio-single-cell-clustering**: PCA、UMAP/tSNE 降维，Leiden/Louvain 聚类
- **scanpy**: 完整的单细胞分析工作流（推荐使用）
- **anndata**: 单细胞数据结构操作

### 3. 细胞注释
- **bio-single-cell-cell-annotation**: 自动细胞类型注释
- **bio-single-cell-markers-annotation**: 基于 marker 基因的注释

### 4. 批次整合
- **bio-single-cell-batch-integration**: Harmony、Scanorama、BBKNN 等批次校正方法

### 5. 下游分析
- **bio-single-cell-trajectory-inference**: 拟时序分析（Monocle、PAGA）
- **bio-single-cell-cell-communication**: 细胞间通讯分析（CellPhoneDB、NicheNet）
- **bio-single-cell-lineage-tracing**: 谱系追踪分析
- **bio-single-cell-perturb-seq**: 扰动实验分析

### 6. 多组学整合
- **bio-single-cell-multimodal-integration**: CITE-seq、多组学数据整合
- **bio-single-cell-scatac-analysis**: scATAC-seq 分析

### 7. 高级工具
- **scvi-tools**: 基于深度学习的单细胞分析
- **scvelo**: RNA 速率分析

### 8. 扰动预测基准测试
- **scgen_meta_benchmark_skill**: 单细胞扰动预测元基准框架，集成 27 种算法（4 大模块）：
  - Module A（生成模型 7 个）: scGen, trVAE, scGAN, CPA, VAE, CVAE, scVI
  - Module B（深度网络 8 个）: MLP, ResNet, Transformer, DenseNet, U-Net, Attention-MLP, Highway, Capsule
  - Module C（经典 ML 6 个）: Random Forest, XGBoost, Ridge, ElasticNet, SVR, KNN
  - Module D（线性方法 6 个）: DE, Mean Shift, Nearest Neighbor, Linear Projection, ComBat, Scanorama
  - 支持 OOD 评估、Optuna 超参优化、R²/MSE/Pearson 等指标
