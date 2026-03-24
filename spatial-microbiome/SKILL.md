---
name: spatial-microbiome
description: >
  空间转录组中微生物与宿主互作分析 skill。基于 Python（Scanpy/Squidpy）生态，
  处理同一空间坐标体系下的宿主 h5ad 和微生物 h5ad 数据，完成共定位分析、
  空间相关性分析、微生物对宿主基因表达的空间影响评估、微生物空间分布模式聚类等任务。
  输出包括可复现的 Python 代码、出版级可视化图表和结构化分析报告。
  当用户提到以下任何关键词时触发此 skill：空间转录组微生物、微生物宿主互作、
  spatial transcriptomics microbiome、host-microbe interaction、微生物 h5ad、
  共定位分析 colocalization、空间相关性、微生物空间分布、
  squidpy 微生物、scanpy 微生物，或用户上传了两个 h5ad 文件并希望做空间层面的
  微生物-宿主关联分析。即使用户只是模糊地提到"分析微生物和宿主的空间关系"，
  也应使用此 skill。
---

# Spatial Microbiome-Host Interaction Analysis

面向生信/计算生物学研究人员的空间转录组微生物-宿主互作分析工具。

## 适用场景

用户拥有同一空间坐标体系下的两个 AnnData（h5ad）对象：
- **宿主 h5ad**：空间转录组数据（基因表达矩阵 + 空间坐标）
- **微生物 h5ad**：微生物丰度矩阵（物种/OTU/ASV × spots + 相同空间坐标）

目标是在空间维度上揭示微生物与宿主细胞/基因表达之间的关系。

## 前置准备

### 微生物分类层级解析

微生物特征名称包含分类层级前缀，必须先解析：

| 前缀 | 层级 | 英文 |
|------|------|------|
| `p__` | 门 | Phylum |
| `c__` | 纲 | Class |
| `o__` | 目 | Order |
| `f__` | 科 | Family |
| `g__` | 属 | Genus |
| `s__` | 种 | Species |

```python
import re

TAXON_LEVELS = {
    'p__': 'Phylum', 'c__': 'Class', 'o__': 'Order',
    'f__': 'Family', 'g__': 'Genus', 's__': 'Species'
}
TAXON_ORDER = ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

def parse_taxon_level(name):
    """从微生物特征名中提取分类层级前缀"""
    for prefix in TAXON_ORDER:
        if name.startswith(prefix):
            return prefix
    return 'unknown'

def split_by_taxon_level(adata_microbe):
    """按分类层级拆分微生物 AnnData，返回 {prefix: adata_subset}"""
    levels = [parse_taxon_level(v) for v in adata_microbe.var_names]
    adata_microbe.var['taxon_level'] = levels
    result = {}
    for prefix in TAXON_ORDER:
        mask = adata_microbe.var['taxon_level'] == prefix
        if mask.sum() > 0:
            result[prefix] = adata_microbe[:, mask].copy()
            print(f"  {TAXON_LEVELS.get(prefix, prefix)}: {mask.sum()} taxa")
    return result
```

后续所有分析（相关性、共定位、差异表达）应按层级分别进行，避免混合不同分类粒度。

### 输入数据校验与空间对齐

宿主和微生物数据来自同一张切片、同一坐标系，但空间分辨率通常不同：
- **宿主数据**：规则网格（如步长 50 的 bin）
- **微生物数据**：更细的不规则坐标（单细胞/亚 bin 级别）

正确的对齐方式是将微生物数据聚合到宿主的网格 bin 中：

```python
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# 1. 加载数据
adata_host = sc.read_h5ad("host.h5ad")
adata_microbe = sc.read_h5ad("microbe.h5ad")

# 2. 校验空间坐标
assert "spatial" in adata_host.obsm, "宿主数据缺少 .obsm['spatial']"
assert "spatial" in adata_microbe.obsm, "微生物数据缺少 .obsm['spatial']"

# 3. 检测宿主网格步长
host_coords = adata_host.obsm['spatial'].astype(float)
hx_unique = np.sort(np.unique(host_coords[:, 0]))
host_step = np.median(np.diff(hx_unique))
print(f"宿主网格步长: {host_step}")

# 4. 将微生物 spots 分配到最近的宿主 bin
micro_coords = adata_microbe.obsm['spatial'].astype(float)
host_tree = cKDTree(host_coords)
dists, nearest_host_idx = host_tree.query(micro_coords, k=1)

# 5. 按宿主 bin 聚合微生物丰度（求和）
micro_X = adata_microbe.X.toarray() if hasattr(adata_microbe.X, 'toarray') else adata_microbe.X
n_host = adata_host.n_obs
n_taxa = adata_microbe.n_vars
aggregated = np.zeros((n_host, n_taxa))
for i in range(len(nearest_host_idx)):
    aggregated[nearest_host_idx[i]] += micro_X[i]

# 6. 创建对齐后的微生物 AnnData（与宿主 spots 一一对应）
adata_micro_aligned = sc.AnnData(
    X=aggregated,
    obs=adata_host.obs.copy(),
    var=adata_microbe.var.copy(),
    obsm={'spatial': host_coords.copy()}
)
print(f"对齐后: {adata_micro_aligned.n_obs} spots × {adata_micro_aligned.n_vars} taxa")
print(f"非零 spots: {(aggregated.sum(axis=1) > 0).sum()}")
```

对齐后，宿主和微生物数据具有完全相同的 spots（行），可以直接做逐 spot 的相关性和差异分析。

### 数据预处理建议

根据数据状态决定是否需要以下步骤（先检查数据是否已处理过，不要重复处理）：

- **宿主数据**：检查是否已做过 normalize_total、log1p、highly_variable_genes 等
- **微生物数据**：检查是否需要 CLR 变换或相对丰度归一化（微生物组合数据的组成性偏差很重要）
- **空间邻域图**：对齐后的数据使用相同的空间坐标构建邻域图

```python
# 如果宿主数据尚未预处理
if "log1p" not in adata_host.uns:
    sc.pp.normalize_total(adata_host, target_sum=1e4)
    sc.pp.log1p(adata_host)

# 微生物数据 CLR 变换（应对组成性数据，在对齐后的数据上做）
from scipy.stats import gmean
def clr_transform(adata):
    """对微生物丰度做 CLR 变换"""
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
    X = X + 0.5  # pseudocount
    gm = gmean(X, axis=1, keepdims=True)
    X_clr = np.log(X / gm)
    adata.layers["clr"] = X_clr
    return adata

adata_micro_aligned = clr_transform(adata_micro_aligned)
```

---

## 核心分析模块

根据用户需求，从以下四个模块中选择合适的组合。如果用户没有明确指定，建议按模块 1→4 的顺序依次进行，因为后续分析常依赖前面的结果。

### 模块 1：微生物空间分布模式与聚类

目标：了解微生物在组织空间上的分布规律。

关键分析步骤：

1. **空间可视化**：在组织坐标上绘制各微生物的丰度分布热图
2. **空间自相关**：用 Moran's I 检验微生物丰度是否有空间聚集趋势
3. **空间聚类**：识别微生物组成相似的空间区域

```python
# 空间可视化 - 绘制 top 微生物的空间分布
top_microbes = adata_microbe.var_names[:6]  # 或按丰度排序取 top
sq.pl.spatial_scatter(
    adata_microbe, color=top_microbes, 
    ncols=3, figsize=(15, 10), cmap="Reds",
    title=[f"{m} 空间分布" for m in top_microbes]
)

# Moran's I 空间自相关
sq.gr.spatial_autocorr(adata_microbe, mode="moran", n_jobs=-1)
morans_results = adata_microbe.uns["moranI"]
sig_microbes = morans_results[morans_results["pval_norm"] < 0.05].sort_values("I", ascending=False)
print(f"具有显著空间聚集的微生物: {len(sig_microbes)}")
print(sig_microbes.head(10))

# 空间聚类 - 识别微生物空间 niche
sc.pp.pca(adata_microbe, n_comps=min(20, adata_microbe.n_vars - 1))
sc.pp.neighbors(adata_microbe, use_rep="X_pca")
sc.tl.leiden(adata_microbe, resolution=0.5, key_added="microbe_niche")
sq.pl.spatial_scatter(adata_microbe, color="microbe_niche", figsize=(8, 8))
```

> 详细的可视化定制和参数调优见 `references/visualization.md`

### 模块 2：共定位分析（微生物-宿主细胞空间邻近）

目标：量化微生物与特定宿主细胞类型在空间上的邻近程度。

前提条件：宿主数据需要有细胞类型注释（存放在 `adata_host.obs` 的某个列中）。如果没有，先帮用户做细胞类型注释或让用户提供。

核心思路：在空间邻域图中，检验某种微生物高丰度区域是否与特定宿主细胞类型显著邻近。

```python
# 将微生物 niche 标签和宿主细胞类型合并到同一个 AnnData
adata_combined = adata_host.copy()
# 添加微生物丰度信息
for microbe in adata_microbe.var_names:
    col = adata_microbe[:, microbe].X.toarray().flatten() \
        if hasattr(adata_microbe.X, "toarray") \
        else adata_microbe[:, microbe].X.flatten()
    adata_combined.obs[f"microbe_{microbe}"] = col

# 基于丰度将 spots 标记为微生物阳性/阴性
def label_positive_spots(adata, microbe_col, threshold_quantile=0.75):
    """基于分位数阈值标记微生物阳性 spots"""
    vals = adata.obs[microbe_col]
    threshold = vals.quantile(threshold_quantile)
    adata.obs[f"{microbe_col}_pos"] = (vals > threshold).astype(str)
    return adata

# Neighborhood enrichment analysis
# 检验微生物阳性区域与宿主细胞类型的空间邻近富集
cell_type_col = "cell_type"  # 根据实际列名调整
sq.gr.spatial_neighbors(adata_combined, coord_type="generic", n_neighs=6)

# 对每种感兴趣的微生物做邻域富集分析
target_microbe = "microbe_Bacteroides"  # 示例
adata_combined = label_positive_spots(adata_combined, target_microbe)

sq.gr.nhood_enrichment(
    adata_combined, 
    cluster_key=f"{target_microbe}_pos"
)
sq.pl.nhood_enrichment(
    adata_combined, 
    cluster_key=f"{target_microbe}_pos",
    title=f"{target_microbe} 与宿主细胞类型共定位"
)
```

对于更精细的共定位分析，考虑使用排列检验：

```python
def permutation_colocalization(adata, microbe_col, celltype_col, 
                                celltype, n_perm=1000):
    """排列检验评估微生物与特定细胞类型的空间共定位显著性"""
    from scipy.sparse import issparse
    
    # 实际观测：微生物阳性 spots 中目标细胞类型的比例
    pos_mask = adata.obs[f"{microbe_col}_pos"] == "True"
    observed = (adata.obs.loc[pos_mask, celltype_col] == celltype).mean()
    
    # 排列分布
    perm_scores = []
    for _ in range(n_perm):
        shuffled = np.random.permutation(adata.obs[celltype_col].values)
        perm_scores.append((shuffled[pos_mask] == celltype).mean())
    
    perm_scores = np.array(perm_scores)
    p_value = (perm_scores >= observed).mean()
    z_score = (observed - perm_scores.mean()) / (perm_scores.std() + 1e-10)
    
    return {"observed": observed, "p_value": p_value, "z_score": z_score}
```

### 模块 3：空间相关性分析（微生物丰度与宿主基因表达）

目标：发现哪些宿主基因的表达在空间上与微生物丰度相关。

核心方法：
- **全局相关性**：Spearman 相关（不受线性假设约束）
- **空间加权相关性**：考虑空间自相关后的校正相关性
- **局部空间相关性**：Lee's L 统计量（双变量空间自相关）

```python
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

def spatial_correlation_screen(adata_host, adata_microbe, 
                                target_microbe, 
                                gene_subset=None,
                                method="spearman"):
    """
    筛选与目标微生物丰度空间相关的宿主基因。
    
    Parameters
    ----------
    adata_host : AnnData - 宿主表达数据
    adata_microbe : AnnData - 微生物丰度数据
    target_microbe : str - 目标微生物名
    gene_subset : list, optional - 只检测这些基因（加速计算）
    method : str - "spearman" 或 "pearson"
    
    Returns
    -------
    pd.DataFrame - 基因名、相关系数、p值、校正p值
    """
    microbe_abundance = adata_microbe[:, target_microbe].X
    if hasattr(microbe_abundance, "toarray"):
        microbe_abundance = microbe_abundance.toarray().flatten()
    else:
        microbe_abundance = microbe_abundance.flatten()
    
    genes = gene_subset if gene_subset else adata_host.var_names
    results = []
    
    for gene in genes:
        gene_expr = adata_host[:, gene].X
        if hasattr(gene_expr, "toarray"):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = gene_expr.flatten()
        
        if method == "spearman":
            corr, pval = spearmanr(microbe_abundance, gene_expr)
        else:
            from scipy.stats import pearsonr
            corr, pval = pearsonr(microbe_abundance, gene_expr)
        
        results.append({"gene": gene, "correlation": corr, "pval": pval})
    
    df = pd.DataFrame(results)
    df["pval_adj"] = multipletests(df["pval"], method="fdr_bh")[1]
    df = df.sort_values("pval_adj")
    
    return df

# 使用示例：筛选与 Bacteroides 空间相关的宿主基因
# 建议先限制到 highly variable genes 以加速
if "highly_variable" in adata_host.var.columns:
    hvg = adata_host.var_names[adata_host.var["highly_variable"]]
else:
    sc.pp.highly_variable_genes(adata_host, n_top_genes=2000)
    hvg = adata_host.var_names[adata_host.var["highly_variable"]]

corr_df = spatial_correlation_screen(
    adata_host, adata_microbe,
    target_microbe="Bacteroides",
    gene_subset=hvg.tolist()
)

sig_genes = corr_df[corr_df["pval_adj"] < 0.05]
print(f"显著相关基因数: {len(sig_genes)}")
```

对显著基因做通路富集分析：

```python
# 正相关基因和负相关基因分开分析
pos_genes = sig_genes[sig_genes["correlation"] > 0]["gene"].tolist()
neg_genes = sig_genes[sig_genes["correlation"] < 0]["gene"].tolist()

# 使用 gseapy 或 gprofiler 做富集（需要安装）
# pip install gseapy
import gseapy as gp

if pos_genes:
    enr_pos = gp.enrichr(
        gene_list=pos_genes,
        gene_sets=["GO_Biological_Process_2021", "KEGG_2021_Human"],
        organism="Human",
        outdir=None
    )
    print("=== 正相关基因富集通路 ===")
    print(enr_pos.results.head(10))
```

### 模块 4：微生物对宿主基因表达的空间影响

目标：在空间上建模微生物对宿主基因表达（尤其是免疫通路）的影响，超越简单相关性。

方法包括：
- **空间分箱对比**：按微生物丰度高/低区域，对比宿主基因表达差异
- **空间回归模型**：考虑空间自相关的回归分析
- **配体-受体空间互作分析**（如果适用）

```python
from scipy.stats import mannwhitneyu

def spatial_impact_analysis(adata_host, adata_microbe, 
                             target_microbe, 
                             gene_set=None,
                             quantile_threshold=0.75):
    """
    按微生物丰度高/低区域对宿主基因做差异表达分析。
    
    类似于空间上的 "伪 bulk" 差异分析：
    - 高丰度区：微生物丰度 > 指定分位数的 spots
    - 低丰度区：微生物丰度 < (1 - 指定分位数) 分位数的 spots
    """
    microbe_vals = adata_microbe[:, target_microbe].X
    if hasattr(microbe_vals, "toarray"):
        microbe_vals = microbe_vals.toarray().flatten()
    else:
        microbe_vals = microbe_vals.flatten()
    
    high_thresh = np.quantile(microbe_vals, quantile_threshold)
    low_thresh = np.quantile(microbe_vals, 1 - quantile_threshold)
    
    high_mask = microbe_vals > high_thresh
    low_mask = microbe_vals < low_thresh
    
    genes = gene_set if gene_set else adata_host.var_names
    results = []
    
    for gene in genes:
        expr = adata_host[:, gene].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        
        high_expr = expr[high_mask]
        low_expr = expr[low_mask]
        
        if len(high_expr) > 5 and len(low_expr) > 5:
            stat, pval = mannwhitneyu(high_expr, low_expr, alternative="two-sided")
            log2fc = np.log2((high_expr.mean() + 1e-6) / (low_expr.mean() + 1e-6))
            results.append({
                "gene": gene, "log2FC": log2fc, "pval": pval,
                "mean_high": high_expr.mean(), "mean_low": low_expr.mean()
            })
    
    df = pd.DataFrame(results)
    df["pval_adj"] = multipletests(df["pval"], method="fdr_bh")[1]
    df = df.sort_values("pval_adj")
    return df

# 聚焦免疫相关基因（示例基因列表，实际应根据物种调整）
immune_genes = ["TNF", "IL1B", "IL6", "IFNG", "IL10", "CXCL8", 
                "TLR2", "TLR4", "NOD2", "MUC2", "DEFA5", "REG3A"]
# 过滤掉数据中不存在的基因
immune_genes = [g for g in immune_genes if g in adata_host.var_names]

impact_df = spatial_impact_analysis(
    adata_host, adata_microbe,
    target_microbe="Bacteroides",
    gene_set=immune_genes
)
print(impact_df)
```

---

## 可视化规范

所有图表应达到出版质量，遵循以下规范：

### 基本设置

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

# 出版级设置
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.figsize": (8, 6),
    "savefig.bbox": "tight",
    "savefig.transparent": True,
})
```

### 必备图表类型

1. **空间分布图**：微生物丰度叠加在组织坐标上（`sq.pl.spatial_scatter`）
2. **相关性热图**：微生物-基因空间相关性矩阵（`sns.clustermap`）
3. **共定位评分图**：邻域富集 z-score 热图
4. **Volcano plot**：微生物高/低丰度区域的差异基因
5. **空间 niche 图**：微生物空间聚类 + 宿主细胞类型叠加

```python
import seaborn as sns

# Volcano plot 示例
def plot_volcano(df, fc_col="log2FC", pval_col="pval_adj", 
                 fc_thresh=0.5, pval_thresh=0.05, title="Volcano Plot"):
    fig, ax = plt.subplots(figsize=(8, 6))
    
    df["-log10p"] = -np.log10(df[pval_col].clip(lower=1e-300))
    
    # 分类着色
    sig_up = (df[fc_col] > fc_thresh) & (df[pval_col] < pval_thresh)
    sig_down = (df[fc_col] < -fc_thresh) & (df[pval_col] < pval_thresh)
    
    ax.scatter(df.loc[~(sig_up | sig_down), fc_col], 
               df.loc[~(sig_up | sig_down), "-log10p"], 
               c="grey", alpha=0.5, s=10, label="Not sig.")
    ax.scatter(df.loc[sig_up, fc_col], df.loc[sig_up, "-log10p"], 
               c="#e74c3c", alpha=0.7, s=15, label="Up")
    ax.scatter(df.loc[sig_down, fc_col], df.loc[sig_down, "-log10p"], 
               c="#3498db", alpha=0.7, s=15, label="Down")
    
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(title)
    ax.legend()
    ax.axhline(-np.log10(pval_thresh), ls="--", c="grey", lw=0.8)
    ax.axvline(fc_thresh, ls="--", c="grey", lw=0.8)
    ax.axvline(-fc_thresh, ls="--", c="grey", lw=0.8)
    
    plt.tight_layout()
    return fig
```

> 更多可视化模板见 `references/visualization.md`

---

## 输出规范

每次分析结束后，生成以下交付物：

### 1. 可复现代码

以完整 Python 脚本形式输出，包含：
- 所有 import 和环境依赖注释
- 数据加载与校验
- 分析流程（带详细注释）
- 图表生成与保存
- 结果表格导出

脚本头部必须包含依赖说明：

```python
"""
Spatial Microbiome-Host Interaction Analysis
=============================================
Dependencies:
    pip install scanpy squidpy anndata numpy pandas scipy 
    pip install matplotlib seaborn statsmodels gseapy
"""
```

### 2. 可视化图表

- 保存为 PNG（300 dpi）和 PDF 两种格式
- 文件名有意义：`colocalization_Bacteroides_immune_cells.png`
- 每张图附简短标题和说明

### 3. 结构化分析报告

报告应包含以下章节（用 Markdown 或 Word 格式输出）：

```
# 空间微生物-宿主互作分析报告

## 1. 数据概览
- 样本信息、spots 数、微生物与基因特征数
- 数据质量评估

## 2. 微生物空间分布
- 关键微生物的空间模式（附图）
- 空间自相关检验结果
- 微生物空间 niche 划分

## 3. 微生物-宿主共定位
- 共定位分析结果（附热图）
- 显著共定位的微生物-细胞类型对

## 4. 空间相关性
- 与微生物丰度显著相关的宿主基因
- 通路富集结果

## 5. 微生物对宿主基因的空间影响
- 微生物高/低丰度区的差异表达基因
- 免疫通路分析
- 关键发现总结

## 6. 方法学说明
- 统计方法、校正方式、软件版本
```

---

## 注意事项

- **组成性数据问题**：微生物丰度是组成性数据，直接做相关分析会产生虚假相关。始终提醒用户使用 CLR 或其他组成性校正方法，并在报告中说明。
- **多重检验校正**：涉及大规模检验时必须做 FDR 校正，默认使用 Benjamini-Hochberg。
- **空间自相关的影响**：空间数据点不独立，简单相关性检验的 p 值可能过于乐观。在报告中注明这一局限，或使用空间校正方法。
- **稀疏性处理**：微生物数据通常非常稀疏，分析前需讨论零值处理策略（pseudocount、零膨胀模型等）。
- **可复现性**：代码中记录随机种子（`np.random.seed(42)`）和软件版本。
