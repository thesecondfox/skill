# Changelog: Local Reference Registry & Offline Enrichment Support
> Date: 2026-04-05
> Author: thesecondfox

## Summary

新增 `local-reference-registry` skill，实现本地参考数据统一管理。修改 4 个下游 skill 支持 local-first 模式，解决 clusterProfiler 联网问题。

## Problem

clusterProfiler 的 `enrichGO()`, `enrichKEGG()`, `gseGO()`, `gseKEGG()` 等函数默认依赖网络连接获取基因集数据（OrgDb, KEGG API）。在无网络环境（如 HPC 集群计算节点）或网络不稳定时，分析无法运行。用户已下载全套 GMT 文件（人和小鼠），但缺乏统一的路径管理机制。

## Solution

引入「参考数据注册中心」模式：

```
~/.bioref_config.yaml      ← 统一配置文件（不绑定任何特定目录）
    ↓ 被读取
┌─────────────────────────────────────────────────┐
│  bio-pathway-analysis-go-enrichment             │
│  bio-pathway-analysis-gsea                      │  → enricher() / GSEA() + TERM2GENE
│  bio-pathway-analysis-kegg-pathways             │     替代 enrichGO() / enrichKEGG()
│  bio-single-cell-cell-communication (CellChat)  │  → 支持本地 RDS 配受体库
└─────────────────────────────────────────────────┘
```

## Changes

### New File

| File | Description |
|------|-------------|
| `Common_Skills/local-reference-registry/SKILL.md` | 参考数据注册中心 skill，包含路径扫描、GMT 自动分类、配置文件生成、R/Python 加载代码 |

### Modified Files

| File | Change |
|------|--------|
| `Omics_Domains/Pathway_Analysis/bio-pathway-analysis-go-enrichment/SKILL.md` | 添加 `depends_on: local-reference-registry`；新增 Local-First Mode 章节，使用 `enricher()` + 本地 GMT 替代 `enrichGO()` |
| `Omics_Domains/Pathway_Analysis/bio-pathway-analysis-gsea/SKILL.md` | 添加 `depends_on: local-reference-registry`；新增 Local-First Mode 章节，使用 `GSEA()` + 本地 GMT 替代 `gseGO()` / `gseKEGG()` |
| `Omics_Domains/Pathway_Analysis/bio-pathway-analysis-kegg-pathways/SKILL.md` | 添加 `depends_on: local-reference-registry`；新增 Local-First Mode 章节，使用 `enricher()` + 本地 KEGG GMT 替代 `enrichKEGG()` |
| `Omics_Domains/Single_Cell/bio-single-cell-cell-communication/SKILL.md` | 添加 `depends_on: local-reference-registry`；新增 Local-First Mode 章节，支持从本地 RDS 加载自定义配受体数据库 |
| `Omics_Domains/Pathway_Analysis/README.md` | 新增"参考数据配置"章节和"离线模式说明" |

## Technical Details

### Configuration File Format

位置：`~/.bioref_config.yaml`（用户 home 目录，不绑定任何工具目录）

```yaml
ref_root: /path/to/your/reference/data
gmt_files:
  human:
    hallmark: h.all.v2024.1.Hs.symbols.gmt
    kegg: c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt
    go_bp: c5.go.bp.v2024.1.Hs.symbols.gmt
    go_mf: c5.go.mf.v2024.1.Hs.symbols.gmt
    go_cc: c5.go.cc.v2024.1.Hs.symbols.gmt
    reactome: c2.cp.reactome.v2024.1.Hs.symbols.gmt
    ...
  mouse:
    hallmark: mh.all.v2024.1.Mm.symbols.gmt
    ...
cellchat_db:
  human: CellChatDB.human.rds
  mouse: CellChatDB.mouse.rds
```

### API Change Pattern

| Before (online) | After (local) |
|-----------------|---------------|
| `enrichGO(gene, OrgDb=org.Hs.eg.db, ont="BP")` | `enricher(gene, TERM2GENE=read.gmt(local_go_bp))` |
| `enrichKEGG(gene, organism="hsa")` | `enricher(gene, TERM2GENE=read.gmt(local_kegg))` |
| `gseGO(geneList, OrgDb=...)` | `GSEA(geneList, TERM2GENE=read.gmt(local_go_bp))` |
| `gseKEGG(geneList, organism=...)` | `GSEA(geneList, TERM2GENE=read.gmt(local_kegg))` |
| `CellChatDB.human` (built-in) | `readRDS(local_cellchat_db)` |

### Behavior Rules

1. 所有下游 skill 在执行富集分析前先检查 `~/.bioref_config.yaml`
2. 配置存在 → 自动使用本地模式（无网络依赖）
3. 配置不存在 → 询问用户是否有本地数据，无则回退到联网模式
4. 路径不硬编码，换机器只需修改配置文件中的 `ref_root`

### Supported GMT Categories

| Key | MSigDB | Description |
|-----|--------|-------------|
| hallmark | H | Hallmark gene sets |
| kegg | C2:KEGG | KEGG pathways |
| reactome | C2:Reactome | Reactome pathways |
| go_bp | C5:GO:BP | GO Biological Process |
| go_mf | C5:GO:MF | GO Molecular Function |
| go_cc | C5:GO:CC | GO Cellular Component |
| wikipathways | C2:WikiPathways | WikiPathways |
| tft | C3:TFT | Transcription factor targets |
| immunesigdb | C7 | Immunologic signatures |
| oncogenic | C6 | Oncogenic signatures |

## Backward Compatibility

- 所有修改向后兼容：无 `~/.bioref_config.yaml` 时行为与修改前完全一致
- 原有的 `enrichGO()` / `enrichKEGG()` 联网代码保留在 "Online Mode" 章节
- 不影响任何已有的分析脚本
