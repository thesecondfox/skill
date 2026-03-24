---
name: cell-annotation
description: Automated cell type annotation for single-cell RNA-seq data using multi-strategy consensus approach with marker-based, reference-based, and LLM-assisted methods
---

# Cell Annotation Skill

## Overview

This skill provides end-to-end automated cell type annotation for single-cell RNA-seq data. It combines multiple annotation strategies (marker-based, reference-based, LLM-assisted) through a rigorous multi-agent deliberation framework to produce high-confidence, biologically validated annotations.

## When to Use This Skill

Trigger this skill when the user mentions:
- "cell annotation"
- "cell type annotation"
- "annotate cells"
- "identify cell types"
- "cell type identification"
- "label clusters"
- "what cell types are these"
- "PBMC annotation" (or any tissue + annotation)
- "marker genes"
- "SingleR" / "CellTypist" / "Azimuth" (reference-based tools)

## Prerequisites

**Required**:
- Single-cell data in h5ad format (AnnData object)
- Data must be preprocessed: normalized, log-transformed, clustered
- Clustering results (e.g., Leiden or Louvain clusters)

**Optional**:
- Tissue type (if unknown, will be inferred)
- Custom marker gene list
- User-specified reference dataset

## Skill Architecture

This skill operates in **5 phases**, each with mandatory stop points for user approval:

```
Phase 0: Data Preparation & QC Check
Phase 1: Annotation Strategy Selection (Multi-Agent Deliberation)
Phase 2: Marker Gene Library Preparation
Phase 3: Multi-Strategy Annotation Execution
Phase 4: Results Integration & Validation
Phase 5: Verification & Export
```

### Phase 0: Data Preparation & QC Check

**Purpose**: Validate input data quality and preprocessing status.

**Checks**:
- [ ] Data format (h5ad)
- [ ] Preprocessing status (normalized, log-transformed)
- [ ] Clustering exists
- [ ] Gene names standardized
- [ ] No critical quality issues

**Output**: QC report, proceed/fix decision

---

### Phase 1: Annotation Strategy Selection

**Purpose**: Multi-agent deliberation to determine optimal annotation approach.

**Agents**:
1. **Tissue Biology Expert** — Validates tissue type, lists expected cell types, checks clustering resolution
2. **Annotation Method Specialist** — Recommends methods (marker-based/reference-based/hybrid), evaluates data suitability
3. **Quality Controller** — Audits data quality, identifies risks, flags potential issues

**Deliberation Protocol**:
- Max 5 rounds of discussion
- All agents participate in every round
- Convergence when all agents approve (PASS) or conditions met (CONDITIONAL)
- If no convergence after 5 rounds → escalate to user

**Output**: `annotation-strategy.yaml` with approved method, expected cell types, quality assessment

**STOP**: Wait for user approval before Phase 2

---

### Phase 2: Marker Gene Library Preparation

**Purpose**: Build comprehensive, tissue-specific marker gene library.

**Steps**:
1. **Database Query** (parallel subagent)
   - CellMarker
   - PanglaoDB
   - Tissue-specific databases (ImmGen, Allen Brain, etc.)

2. **User Marker Integration**
   - Load custom markers (if provided)
   - Validate gene names
   - Merge with database markers

3. **Detection Validation** (subagent)
   - Check which markers exist in dataset
   - Resolve synonyms (via MyGene.info)
   - Calculate detection rates

4. **Prioritization**
   - Classify: canonical vs auxiliary
   - Score by specificity
   - Rank by reliability

5. **Conflict Resolution**
   - Identify ambiguous markers (expressed in multiple types)
   - Create marker combinations for disambiguation
   - Document logic

**Output**: `marker_gene_library.yaml`, detection report

**Quality Gate**:
- Average detection rate ≥70%
- No cell type with <50% detection (or user acknowledges risk)

**STOP**: Wait for user approval before Phase 3

---

### Phase 3: Multi-Strategy Annotation Execution

**Purpose**: Run multiple annotation methods in parallel, integrate via consensus voting.

**Methods**:

#### Method A: Marker-Based Annotation (always runs)
- Compute cluster-level marker scores
- 3 scoring algorithms: weighted / binary / rank
- Assign cell type with highest score
- Confidence = gap between top 2 scores

#### Method B: Reference-Based Annotation (if strategy includes)
- Tools: SingleR / CellTypist / Azimuth
- Transfer labels from reference atlas
- Majority voting at cluster level
- Cross-validate with marker expression

#### Method C: LLM-Assisted Annotation (for ambiguous clusters)
- Use GPT-4/Claude to interpret top marker genes
- Provide biological reasoning
- Confidence: HIGH/MEDIUM/LOW

**Consensus Integration**:
- Weighted voting (confidence-weighted)
- Conflict resolution rules:
  1. Marker-based high confidence (≥0.8) → trust it
  2. Check actual marker expression for candidates
  3. Flag for manual review if unresolved

**Confidence Calibration**:
- FULL_AGREEMENT → boost 10%
- STRONG_MAJORITY → no change
- WEAK_MAJORITY → reduce 10%
- DISAGREEMENT → reduce 30%

**Output**: `consensus_annotations.yaml`, manual review list

**Quality Gate**:
- All clusters annotated
- ≥80% clusters with confidence >0.7
- Manual review list reasonable (<20%)

**STOP**: Wait for user to review flagged clusters

---

### Phase 4: Results Integration & Validation

**Purpose**: Apply annotations, validate, generate figures.

**Steps**:
1. **Apply to Data**
   - Add `cell_type` column to adata.obs
   - Add `annotation_confidence`
   - Add `agreement_status`

2. **Final Validation**
   - Cell type proportion check (vs expected from Phase 1)
   - Marker expression validation
   - Cross-cluster consistency

3. **Visualizations** (publication-quality)
   - UMAP with cell types
   - UMAP with confidence
   - Cell type proportions
   - Marker heatmap
   - Confidence distribution

4. **Evidence Report**
   - Per-cluster annotation evidence
   - Method comparison
   - Quality metrics

5. **Export**
   - Annotated h5ad
   - CSV tables
   - Figures (PDF, 300 DPI)

**Output**: Annotated data, figures, evidence report

**STOP**: Wait for user approval before Phase 5

---

### Phase 5: Verification & Export

**Purpose**: Generate reproducible pipeline, documentation, and final deliverables.

**Steps**:
1. **Reproducibility Package**
   - Standalone Python script
   - Environment documentation (requirements.txt, environment.yaml)
   - README with instructions

2. **Quality Metrics**
   - Annotation coverage
   - Confidence distribution
   - Marker consistency
   - Cluster purity
   - **Quality Grade**: EXCELLENT / GOOD / ACCEPTABLE / NEEDS_REVIEW

3. **Downstream Exports**
   - Seurat (RDS)
   - CellxGene (h5ad)
   - GEO submission package

4. **Documentation**
   - Methods section text (for paper)
   - Figure legends
   - Supplementary tables

5. **Final Package**
   - Archive all results (tar.gz)

**Output**: Complete deliverables package

---

## Domain Sanity Check (Enforced Throughout)

<IRON-LAW>
Before reporting ANY annotation results, perform domain sanity check:

**Check 1: Cell type proportions**
- Example: PBMC should have ~60-70% T cells, not 5%
- Flag if observed proportion outside expected range

**Check 2: Marker expression logic**
- Example: CD3+ but CD19+ in same cluster → likely doublet
- Flag if canonical markers not expressed

**Check 3: Rare cell types**
- Clusters <50 cells annotated as rare types → verify carefully
- Require higher confidence threshold

**Check 4: Biological plausibility**
- Example: Neurons in blood sample → impossible
- Flag unexpected cell types for tissue

If sanity check fails → debug before reporting. Do NOT report suspicious results as "findings."
</IRON-LAW>

---

## Key Design Principles

### 1. Multi-Agent Deliberation
Every major decision (strategy, method, validation) involves 3 specialized agents with distinct perspectives. This catches blind spots and ensures robustness.

### 2. Evidence-Driven
Every annotation has a complete evidence chain:
- Which markers support it
- Which methods agree
- What the confidence score means
- Why conflicts were resolved a certain way

### 3. Biological Plausibility First
Technical correctness ≠ biological correctness. All results validated against domain knowledge.

### 4. Reproducibility by Default
Every step scripted, environment documented, random seeds fixed. Third parties can reproduce exactly.

### 5. Progressive Refinement
Phase 1 → Phase 2 → Phase 3 is a refinement loop. Early phases inform later phases. Phase 4a exploration can trigger returns to Phase 2/3.

---

## Output Files

```
docs/
├── 01_intake/
│   └── research-anchor.yaml
├── 02_literature/
│   └── paper-list.md
├── 03_plan/
│   ├── annotation-strategy.yaml
│   └── marker_gene_library.yaml
├── 04_annotation/
│   ├── marker_based_results.yaml
│   ├── reference_based_results.yaml
│   └── consensus_annotations.yaml
├── 05_figures/
│   ├── umap_cell_types.pdf
│   ├── umap_confidence.pdf
│   ├── cell_type_proportions.pdf
│   ├── marker_heatmap.pdf
│   └── confidence_distribution.pdf
└── 06_results/
    ├── annotated_data.h5ad
    ├── cell_annotations.csv
    ├── evidence_report.md
    ├── quality_certificate.yaml
    ├── methods_section.txt
    └── README.md
```

---

## Usage Example

```python
# User provides data
adata = sc.read_h5ad("pbmc_data.h5ad")

# Trigger skill
# User: "Please annotate the cell types in this PBMC data"

# Skill executes Phase 0-5 with user approval at each gate
# Final output: annotated_data.h5ad with cell_type labels
```

---

## Quality Guarantees

- **Annotation Coverage**: ≥95% of cells annotated
- **High Confidence**: ≥80% of cells with confidence ≥0.7
- **Biological Validity**: All annotations pass domain sanity check
- **Reproducibility**: Complete pipeline script + environment
- **Documentation**: Methods text ready for publication

---

## Integration with Amplify Framework

This skill follows Amplify's core principles:

1. **Phase-by-phase execution** with mandatory stops
2. **Multi-agent deliberation** for critical decisions
3. **Evidence-based claims** with verification
4. **Domain sanity checks** enforced
5. **Reproducibility-driven** from start

---

## Limitations

- Requires preprocessed, clustered data (does not perform clustering)
- Reference-based methods require suitable reference datasets
- LLM-assisted annotation requires API access (Claude/GPT-4)
- Novel cell types not in databases may need manual curation
- Doublet detection should be performed before annotation

---

## Citation

If you use this skill, please cite:
- Amplify Cell Annotation Skill (2026)
- CellMarker: http://xteam.xbio.top/CellMarker/
- PanglaoDB: https://panglaodb.se/

---

## Support

For issues or questions:
- Check the evidence report: `docs/06_results/evidence_report.md`
- Review the quality certificate: `docs/06_results/quality_certificate.yaml`
- Consult the reproducible pipeline: `code/annotation_pipeline.py`

---

**Version**: 1.0
**Last Updated**: 2026-03-05
**Maintainer**: Amplify Framework
