# Skills 整合报告
> 日期: 2026-03-26
> 提交: 072e714
> 仓库: thesecondfox/skill

## 改动摘要

本次将 2 个顶层散落的 skill 整合进 Two-Tier 目录结构（Common_Skills / Omics_Domains）。

## 具体改动

### 1. bio-aging-clocks → Omics_Domains/Epigenomics/

| 项目 | 内容 |
|------|------|
| 原路径 | `/bio-aging-clocks/SKILL.md` |
| 新路径 | `/Omics_Domains/Epigenomics/bio-aging-clocks/SKILL.md` |
| 归类依据 | 该 skill 基于 DNA 甲基化数据进行表观遗传衰老时钟分析（Horvath, GrimAge, DunedinPACE 等），属于表观基因组学范畴 |
| 上下游关系 | 与 `bio-methylation-analysis-*` 系列形成完整链路：原始甲基化数据处理 → DMR 检测 → 衰老时钟预测 |
| README 更新 | 在 Epigenomics/README.md 第 3 节"DNA 甲基化分析"中新增条目 |

### 2. scgen_meta_benchmark_skill → Omics_Domains/Single_Cell/

| 项目 | 内容 |
|------|------|
| 原路径 | `/scgen_meta_benchmark_skill/` (含 core/, sub_skills/, examples/, data/, results/) |
| 新路径 | `/Omics_Domains/Single_Cell/scgen_meta_benchmark_skill/` |
| 归类依据 | 该 skill 是单细胞扰动预测元基准框架（27 种算法），属于单细胞分析范畴 |
| 上下游关系 | 与 `bio-single-cell-perturb-seq` 互补：perturb-seq 处理扰动实验原始数据，scgen 负责扰动响应预测与多模型对比 |
| README 更新 | 在 Single_Cell/README.md 新增第 8 节"扰动预测基准测试"，列出 4 大模块 27 个模型 |

## 整合原则

1. **按组学归类**: 根据 skill 的核心数据类型和分析目标，归入对应的 Omics_Domains 子目录
2. **Pipeline 串联**: 新 skill 在 README 中按分析流程的逻辑位置插入，而非简单追加
3. **保持结构一致**: 整合后顶层只保留 `Common_Skills/`、`Omics_Domains/`、`README.md`

## 整合后目录状态

```
skills_bioinformatics_only/
├── Common_Skills/              248 个通用工具
├── Omics_Domains/
│   ├── Epigenomics/            26 个 skill（+1 bio-aging-clocks）
│   ├── Single_Cell/            16 个 skill（+1 scgen_meta_benchmark_skill）
│   └── ... 其他 15 个领域
└── README.md
```

## 文件变更清单

- 修改: `Omics_Domains/Epigenomics/README.md`
- 修改: `Omics_Domains/Single_Cell/README.md`
- 重命名: `bio-aging-clocks/` → `Omics_Domains/Epigenomics/bio-aging-clocks/`
- 重命名: `scgen_meta_benchmark_skill/` → `Omics_Domains/Single_Cell/scgen_meta_benchmark_skill/`
- 共计 19 个文件变更，9 行新增
