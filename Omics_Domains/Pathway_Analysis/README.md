# Pathway Analysis & Network Biology

## 分析流程

### 0. 参考数据配置（推荐优先执行）
- **local-reference-registry**（位于 Common_Skills/）: 配置本地 GMT 基因集路径，消除联网依赖。配置后下列富集分析 skill 将自动使用本地数据。

### 1. 富集分析
- **bio-pathway-analysis-go-enrichment**: GO 富集分析（支持本地 GMT 离线模式）
- **bio-pathway-analysis-gsea**: GSEA 基因集富集（支持本地 GMT 离线模式）
- **bio-pathway-analysis-kegg-pathways**: KEGG 通路分析（支持本地 GMT 离线模式）
- **bio-pathway-analysis-reactome-pathways**: Reactome 通路分析
- **bio-pathway-analysis-wikipathways**: WikiPathways 分析
- **bio-pathway-analysis-enrichment-visualization**: 富集结果可视化

### 2. 基因调控网络
- **bio-gene-regulatory-networks-coexpression-networks**: 共表达网络
- **bio-gene-regulatory-networks-scenic-regulons**: SCENIC 调控子分析

### 3. 工具库
- **networkx**: 网络分析库

### 离线模式说明

当网络不可用或需要使用自定义基因集时，请先配置 `local-reference-registry`:
1. 准备本地 GMT 文件（MSigDB 格式，人/小鼠）
2. 运行 `local-reference-registry` 的初始化脚本扫描并注册路径
3. 生成 `~/.bioref_config.yaml` 配置文件
4. 之后 GO/KEGG/GSEA 等 skill 会自动使用本地数据
