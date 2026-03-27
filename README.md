<p align="center">
  <h1 align="center">Bioinformatics Skills Collection</h1>
  <p align="center">
    A curated library of 492 AI-agent skills for computational biology and bioinformatics.
    <br />
    <a href="#omics-domains"><strong>Explore Domains</strong></a>
    &middot;
    <a href="#quick-start"><strong>Quick Start</strong></a>
    &middot;
    <a href="AI_AGENT_GUIDE.md"><strong>AI Agent Guide</strong></a>
  </p>
</p>

<p align="center">
  <a href="https://github.com/thesecondfox/skill/stargazers"><img src="https://img.shields.io/github/stars/thesecondfox/skill?style=flat-square&color=yellow" alt="Stars"></a>
  <a href="https://github.com/thesecondfox/skill/network/members"><img src="https://img.shields.io/github/forks/thesecondfox/skill?style=flat-square" alt="Forks"></a>
  <a href="https://github.com/thesecondfox/skill/issues"><img src="https://img.shields.io/github/issues/thesecondfox/skill?style=flat-square" alt="Issues"></a>
  <a href="LICENSE"><img src="https://img.shields.io/github/license/thesecondfox/skill?style=flat-square" alt="License"></a>
  <a href="https://github.com/thesecondfox/skill/commits/main"><img src="https://img.shields.io/github/last-commit/thesecondfox/skill?style=flat-square" alt="Last Commit"></a>
  <img src="https://img.shields.io/badge/skills-492-blue?style=flat-square" alt="Skills Count">
  <img src="https://img.shields.io/badge/omics_domains-17-green?style=flat-square" alt="Domains">
</p>

---

## Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Omics Domains](#omics-domains)
- [Common Skills](#common-skills)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [For AI Agents](#for-ai-agents)
- [Contributing](#contributing)
- [Star History](#star-history)
- [License](#license)

## Overview

This repository provides **492 curated bioinformatics skills** designed for AI coding agents (Claude Code, Cursor, etc.). Each skill is a structured instruction set that teaches an AI agent how to perform a specific bioinformatics task — from raw data processing to publication-quality visualization.

### Key Features

- **17 omics domains** covering the full spectrum of modern bioinformatics
- **Pipeline-oriented** README in each domain, organized by analysis lifecycle
- **Two-Tier architecture** separating cross-domain tools from domain-specific analyses
- **30+ database connectors** for PubMed, UniProt, KEGG, ClinVar, gnomAD and more
- **40+ pre-built workflows** for RNA-seq, ATAC-seq, GWAS, metagenomics and more
- **AI Agent Guide** with routing table and orchestration protocol

## Architecture

```
skill/
├── Common_Skills/          248 cross-domain tools
│   ├── Database access         (PubMed, UniProt, KEGG, GEO, ...)
│   ├── Visualization           (Circos, heatmaps, volcano, ...)
│   ├── Sequence operations     (I/O, format conversion, intervals)
│   ├── Statistical analysis    (experimental design, multiple testing)
│   ├── Workflow management     (Snakemake, Nextflow, CWL, WDL)
│   ├── Scientific writing      (LaTeX, PPTX, Jupyter, Quarto)
│   └── Core libraries          (scanpy, biopython, pysam, rdkit, ...)
│
├── Omics_Domains/          244 domain-specific skills
│   ├── Genomics/               (39 skills)
│   ├── Epigenomics/            (26 skills)
│   ├── Transcriptomics/        (24 skills)
│   ├── ...                     (14 more domains)
│   └── Long_Read/              (4 skills)
│
├── AI_AGENT_GUIDE.md       Instruction for AI agents
├── CHANGELOG_20260326.md   Change log
└── README.md               This file
```

## Omics Domains

| Domain | Skills | Coverage |
|--------|--------|----------|
| **Genomics** | 39 | Variant calling, genome assembly & annotation, population genetics, structural variants |
| **Epigenomics** | 26 | ChIP-seq, ATAC-seq, DNA methylation, Hi-C, aging clocks |
| **Transcriptomics** | 24 | RNA-seq, differential expression, alternative splicing, Ribo-seq, small RNA |
| **Immunology** | 17 | TCR/BCR analysis, flow cytometry, immunoinformatics, neoantigen prediction |
| **Single Cell** | 16 | Preprocessing, clustering, annotation, trajectory, perturbation prediction (27-model benchmark) |
| **Spatial** | 15 | Spatial transcriptomics, spatial proteomics, domain detection, deconvolution |
| **Alignment** | 14 | BWA/STAR/HISAT2 alignment, SAM/BAM operations, MSA |
| **Functional Genomics** | 14 | CRISPR screens, CLIP-seq, Perturb-seq |
| **Metagenomics** | 14 | Taxonomic classification, functional annotation, AMR detection, strain tracking |
| **Clinical** | 10 | Clinical variant interpretation, liquid biopsy, pharmacogenomics, TMB |
| **Pathway Analysis** | 10 | GO/KEGG/Reactome enrichment, GSEA, WikiPathways |
| **Proteomics** | 10 | DIA analysis, protein quantification, PTM, differential abundance |
| **Chemoinformatics** | 9 | Molecular descriptors, ADMET prediction, virtual screening |
| **Machine Learning** | 8 | Omics classifiers, biomarker discovery, survival analysis |
| **Metabolomics** | 8 | XCMS/MS-DIAL preprocessing, metabolite annotation, lipidomics |
| **Structural Biology** | 6 | AlphaFold prediction, PDB operations, molecular dynamics |
| **Long Read** | 4 | Nanopore/HiFi alignment, Clair3 variant calling, methylation detection |

> Each domain directory contains a `README.md` that maps skills to pipeline stages.

## Common Skills

| Category | Examples | Count |
|----------|----------|-------|
| Database Access | PubMed, UniProt, Ensembl, KEGG, ClinVar, gnomAD, GEO, DrugBank | 30+ |
| Visualization | Circos, genome browser, heatmaps, volcano, UpSet, network, interactive | 15+ |
| Sequence Operations | FASTA/FASTQ I/O, format conversion, BED/GTF intervals, expression matrices | 15+ |
| Workflow Management | Snakemake, Nextflow, CWL, WDL + 40 pre-built pipelines | 40+ |
| Core Libraries | scanpy, anndata, biopython, pysam, rdkit, scikit-learn, matplotlib | 19 |
| Scientific Writing | LaTeX posters, PPTX slides, Jupyter/Quarto/RMarkdown reports | 10+ |
| Experimental Design | Power analysis, sample size, batch design, multiple testing correction | 4 |
| Specialized Topics | Epitranscriptomics, ecological genomics, CRISPR design, systems biology | 25+ |

## Quick Start

### 1. Find the right skill

```
User task: "Analyze my scRNA-seq data"

Step 1: Identify domain        → Single_Cell
Step 2: Read pipeline README   → Omics_Domains/Single_Cell/README.md
Step 3: Pick pipeline stage    → preprocessing → clustering → annotation
Step 4: Load skill             → bio-single-cell-preprocessing/SKILL.md
```

### 2. Skill file structure

Each skill contains at minimum a `SKILL.md`:

```
bio-single-cell-clustering/
└── SKILL.md          # Capabilities, workflow, parameters, examples

# Some skills include additional resources:
bio-aging-clocks/
├── SKILL.md          # Core documentation
├── reference.md      # API reference
└── examples/         # Runnable scripts
```

## Installation

Clone the repository:

```bash
git clone git@github.com:thesecondfox/skill.git
```

### For Claude Code

```bash
# Copy to Claude skills directory
cp -r skill/ ~/.claude/skills/
```

### For Cursor

```bash
# Copy to Cursor rules directory
cp -r skill/ ~/.cursor/rules/
```

## For AI Agents

See [AI_AGENT_GUIDE.md](AI_AGENT_GUIDE.md) for:
- Skill resolution protocol and domain routing table
- Loading priority rules
- Multi-skill orchestration patterns
- Common Skills category index

## Contributing

Contributions are welcome! To add a new skill:

1. Fork this repository
2. Create a new directory under the appropriate location:
   - `Common_Skills/` for cross-domain tools
   - `Omics_Domains/<domain>/` for domain-specific skills
3. Add a `SKILL.md` following the existing format
4. Update the domain `README.md` to include your skill in the pipeline
5. Submit a Pull Request

### Skill format

```markdown
# Skill Name

## Overview
Brief description of what this skill does.

## When to Use This Skill
Conditions that trigger this skill.

## Workflow
Step-by-step analysis procedure.

## Parameters
Key parameters and their defaults.

## Examples
Runnable code examples.
```

## Star History

<a href="https://star-history.com/#thesecondfox/skill&Date">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/svg?repos=thesecondfox/skill&type=Date&theme=dark" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/svg?repos=thesecondfox/skill&type=Date" />
   <img alt="Star History Chart" src="https://api.star-history.com/svg?repos=thesecondfox/skill&type=Date" />
 </picture>
</a>

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  Made with Claude Code
</p>
