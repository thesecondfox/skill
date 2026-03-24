# 可视化参考文档

本文档包含空间微生物-宿主互作分析中常用的高级可视化模板。

## 目录

1. 双层空间叠加图
2. 相关性矩阵聚类热图
3. 空间 niche 组成堆叠图
4. 交互式空间图（可选）
5. 配色方案

---

## 1. 双层空间叠加图

同时展示微生物丰度和宿主细胞类型在组织空间上的分布。

```python
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

def plot_dual_spatial(adata_host, adata_microbe, 
                      microbe_name, celltype_col,
                      figsize=(16, 6)):
    """
    左图：宿主细胞类型空间分布
    中图：微生物丰度空间分布
    右图：叠加（细胞类型底图 + 微生物等高线）
    """
    coords = adata_host.obsm["spatial"]
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # 左图：宿主细胞类型
    celltypes = adata_host.obs[celltype_col].astype("category")
    for i, ct in enumerate(celltypes.cat.categories):
        mask = celltypes == ct
        axes[0].scatter(coords[mask, 0], coords[mask, 1], 
                       s=5, alpha=0.6, label=ct)
    axes[0].set_title("宿主细胞类型分布")
    axes[0].legend(fontsize=7, markerscale=2, 
                   bbox_to_anchor=(1, 1), loc="upper left")
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    
    # 中图：微生物丰度
    microbe_vals = adata_microbe[:, microbe_name].X
    if hasattr(microbe_vals, "toarray"):
        microbe_vals = microbe_vals.toarray().flatten()
    else:
        microbe_vals = microbe_vals.flatten()
    
    sc_mid = axes[1].scatter(coords[:, 0], coords[:, 1], 
                              c=microbe_vals, cmap="Reds", 
                              s=5, alpha=0.7)
    plt.colorbar(sc_mid, ax=axes[1], shrink=0.7, label="丰度")
    axes[1].set_title(f"{microbe_name} 空间丰度")
    axes[1].set_aspect("equal")
    axes[1].invert_yaxis()
    
    # 右图：叠加
    for i, ct in enumerate(celltypes.cat.categories):
        mask = celltypes == ct
        axes[2].scatter(coords[mask, 0], coords[mask, 1], 
                       s=3, alpha=0.3)
    # 微生物丰度等高线
    from scipy.interpolate import griddata
    xi = np.linspace(coords[:, 0].min(), coords[:, 0].max(), 200)
    yi = np.linspace(coords[:, 1].min(), coords[:, 1].max(), 200)
    zi = griddata(coords, microbe_vals, (xi[None, :], yi[:, None]), 
                  method="cubic", fill_value=0)
    axes[2].contour(xi, yi, zi, levels=8, cmap="Reds", alpha=0.8, linewidths=1)
    axes[2].set_title(f"叠加: 细胞类型 + {microbe_name}")
    axes[2].set_aspect("equal")
    axes[2].invert_yaxis()
    
    for ax in axes:
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
    
    plt.tight_layout()
    return fig
```

## 2. 相关性矩阵聚类热图

展示多种微生物与多个宿主基因之间的空间相关性全景。

```python
import seaborn as sns

def plot_correlation_clustermap(corr_matrix, pval_matrix=None,
                                 figsize=(12, 10),
                                 pval_thresh=0.05):
    """
    Parameters
    ----------
    corr_matrix : pd.DataFrame
        行=微生物, 列=宿主基因, 值=相关系数
    pval_matrix : pd.DataFrame, optional
        同维度的 p 值矩阵，不显著的格子用 'x' 标记
    """
    # 标注不显著的位置
    annot = None
    if pval_matrix is not None:
        annot = pval_matrix.applymap(
            lambda x: "" if x < pval_thresh else "×"
        )
    
    g = sns.clustermap(
        corr_matrix,
        cmap="RdBu_r", center=0, vmin=-1, vmax=1,
        figsize=figsize,
        annot=annot if annot is not None else False,
        fmt="",
        linewidths=0.5,
        dendrogram_ratio=(0.15, 0.15),
        cbar_kws={"label": "Spearman ρ", "shrink": 0.5}
    )
    g.ax_heatmap.set_xlabel("宿主基因")
    g.ax_heatmap.set_ylabel("微生物")
    g.fig.suptitle("微生物-宿主基因空间相关性", y=1.02, fontsize=14)
    
    return g

# 构建相关性矩阵的辅助函数
def build_correlation_matrix(adata_host, adata_microbe,
                              microbes=None, genes=None):
    """批量计算微生物-基因相关性矩阵"""
    from scipy.stats import spearmanr
    
    if microbes is None:
        microbes = adata_microbe.var_names[:20]
    if genes is None:
        genes = adata_host.var_names[:50]
    
    corr_mat = pd.DataFrame(index=microbes, columns=genes, dtype=float)
    pval_mat = pd.DataFrame(index=microbes, columns=genes, dtype=float)
    
    for m in microbes:
        m_vals = adata_microbe[:, m].X
        if hasattr(m_vals, "toarray"):
            m_vals = m_vals.toarray().flatten()
        else:
            m_vals = m_vals.flatten()
        
        for g in genes:
            g_vals = adata_host[:, g].X
            if hasattr(g_vals, "toarray"):
                g_vals = g_vals.toarray().flatten()
            else:
                g_vals = g_vals.flatten()
            
            corr, pval = spearmanr(m_vals, g_vals)
            corr_mat.loc[m, g] = corr
            pval_mat.loc[m, g] = pval
    
    return corr_mat.astype(float), pval_mat.astype(float)
```

## 3. 空间 niche 组成堆叠图

展示每个微生物空间 niche 内宿主细胞类型的组成差异。

```python
def plot_niche_composition(adata_host, adata_microbe, 
                            niche_key="microbe_niche",
                            celltype_col="cell_type",
                            figsize=(10, 6)):
    """堆叠柱状图：各微生物 niche 中宿主细胞类型比例"""
    
    # 合并信息
    df = pd.DataFrame({
        "niche": adata_microbe.obs[niche_key].values,
        "celltype": adata_host.obs[celltype_col].values
    })
    
    # 计算比例
    comp = df.groupby("niche")["celltype"].value_counts(normalize=True)
    comp = comp.unstack(fill_value=0)
    
    ax = comp.plot.bar(stacked=True, figsize=figsize, 
                        colormap="Set3", edgecolor="white", linewidth=0.5)
    ax.set_ylabel("细胞类型比例")
    ax.set_xlabel("微生物空间 Niche")
    ax.set_title("各微生物 Niche 中宿主细胞类型组成")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    return ax.figure
```

## 4. 配色方案建议

空间微生物分析推荐的配色方案：

| 用途 | 推荐 colormap | 说明 |
|------|--------------|------|
| 微生物丰度 | `Reds`, `YlOrRd` | 从浅到深，直觉上表示"越多越深" |
| 相关性 | `RdBu_r` | 蓝=负相关，红=正相关，白=无相关 |
| 细胞类型 | `Set3`, `tab20` | 离散颜色，区分度高 |
| 差异表达 | `coolwarm` | Volcano plot 中的上/下调 |
| 空间聚类 | `Paired`, `Set2` | 离散颜色，柔和不刺眼 |
| 显著性 | `Greys` | 灰度表示 -log10(p) |
