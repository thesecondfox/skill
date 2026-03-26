# AI Agent Instruction: Bioinformatics Skills Repository

## Purpose
This repository is a structured skill library for AI agents performing bioinformatics tasks. Each skill is a self-contained unit with a `SKILL.md` file describing capabilities, usage patterns, and constraints. When a user requests a bioinformatics analysis, the agent should locate and load the relevant skill(s) from this repository to guide execution.

## Repository Structure

```
├── Common_Skills/       248 cross-domain skills
└── Omics_Domains/       244 domain-specific skills across 17 omics fields
```

## Skill Resolution Protocol

When receiving a bioinformatics task, follow this decision tree:

1. Identify the omics domain (single-cell, genomics, transcriptomics, etc.)
2. Check `Omics_Domains/<domain>/README.md` for the pipeline stage mapping
3. Load the SKILL.md of the matched skill
4. If the task requires cross-domain tools (databases, visualization, statistics), check `Common_Skills/`

## Domain Routing Table

| User intent keywords | Route to |
|---------------------|----------|
| scRNA-seq, single cell, clustering, UMAP, cell type | Omics_Domains/Single_Cell/ |
| perturbation, drug response, scgen, OOD prediction | Omics_Domains/Single_Cell/scgen_meta_benchmark_skill/ |
| RNA-seq, differential expression, DESeq2, edgeR | Omics_Domains/Transcriptomics/ |
| variant calling, GATK, VCF, SNP, indel | Omics_Domains/Genomics/ |
| ChIP-seq, ATAC-seq, methylation, Hi-C | Omics_Domains/Epigenomics/ |
| aging clock, biological age, Horvath, GrimAge | Omics_Domains/Epigenomics/bio-aging-clocks/ |
| spatial transcriptomics, Visium, MERFISH | Omics_Domains/Spatial/ |
| metagenomics, 16S, microbiome, Kraken | Omics_Domains/Metagenomics/ |
| proteomics, mass spec, DIA, TMT | Omics_Domains/Proteomics/ |
| metabolomics, XCMS, lipidomics | Omics_Domains/Metabolomics/ |
| TCR, BCR, immune repertoire, flow cytometry | Omics_Domains/Immunology/ |
| CRISPR screen, MAGeCK, CLIP-seq | Omics_Domains/Functional_Genomics/ |
| pathway, GO, KEGG, GSEA, enrichment | Omics_Domains/Pathway_Analysis/ |
| clinical variant, ClinVar, pharmacogenomics | Omics_Domains/Clinical/ |
| protein structure, AlphaFold, PDB, MD | Omics_Domains/Structural_Biology/ |
| molecular descriptor, ADMET, virtual screening | Omics_Domains/Chemoinformatics/ |
| nanopore, PacBio, long read, HiFi | Omics_Domains/Long_Read/ |
| alignment, BAM, SAM, BWA, STAR | Omics_Domains/Alignment/ |
| biomarker, classifier, survival, ML for omics | Omics_Domains/Machine_Learning/ |

## Skill File Format

Each skill directory contains at minimum:
- `SKILL.md`: Core document describing persona, capabilities, workflow, constraints, and usage examples

Some skills additionally contain:
- `reference.md`: API documentation, function signatures, parameter details
- `examples/`: Runnable scripts demonstrating usage
- `core/`, `sub_skills/`: Modular code for complex skills (e.g., scgen_meta_benchmark_skill)

## Loading Priority Rules

1. Domain-specific skill > Common skill (prefer specialized over general)
2. Pipeline-aware skill > standalone tool (prefer skills that understand upstream/downstream context)
3. Comprehensive library skill > granular operation skill (e.g., prefer `scanpy` over individual preprocessing skills when running a full workflow)
4. If multiple skills apply, load the one matching the user's current pipeline stage (see domain README.md)

## Common Skills Categories

When the task requires cross-domain capabilities, search Common_Skills by category:

| Category | Prefix/Keywords | Count |
|----------|----------------|-------|
| Database access | bio-database-access-*, *-database | 30+ |
| Visualization | bio-data-visualization-*, matplotlib, seaborn, plotly | 15+ |
| Sequence operations | bio-sequence-io-*, bio-sequence-manipulation-* | 15+ |
| Genomic intervals | bio-genome-intervals-* | 7 |
| Expression matrix | bio-expression-matrix-* | 4 |
| Experimental design | bio-experimental-design-* | 4 |
| Workflow management | bio-workflow-management-*, bio-workflows-* | 40+ |
| Reporting | bio-reporting-*, scientific-writing, pdf, pptx | 10+ |
| Ecological genomics | bio-ecological-genomics-* | 6 |
| Epidemiological | bio-epidemiological-genomics-* | 5 |
| Epitranscriptomics | bio-epitranscriptomics-* | 5 |
| Genome engineering | bio-genome-engineering-* | 5 |
| Comparative genomics | bio-comparative-genomics-* | 5 |
| Systems biology | bio-systems-biology-* | 5 |

## Multi-Skill Orchestration

For complex analyses spanning multiple domains:

1. Load the primary domain skill for the core analysis
2. Load Common_Skills for data I/O, visualization, and statistical testing
3. Load secondary domain skills for integrative analyses
4. Follow the pipeline order defined in each domain's README.md

Example: Single-cell RNA-seq with pathway analysis
```
Load: Single_Cell/bio-single-cell-preprocessing   (step 1: QC)
Load: Single_Cell/bio-single-cell-clustering       (step 2: clustering)
Load: Single_Cell/bio-single-cell-cell-annotation  (step 3: annotation)
Load: Pathway_Analysis/bio-pathway-analysis-gsea   (step 4: downstream)
Load: Common_Skills/bio-data-visualization-*       (visualization)
```

## Constraints

- Always read the full SKILL.md before executing any analysis
- Respect hardware requirements stated in skills (e.g., GPU for deep learning models)
- Follow the data format requirements (e.g., AnnData for single-cell, VCF for variants)
- When a skill specifies version constraints for dependencies, honor them
- Do not mix incompatible tools from different skills without checking compatibility
