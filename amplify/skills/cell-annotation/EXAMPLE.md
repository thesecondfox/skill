# Cell Annotation Skill - Complete Example

## Scenario
User has PBMC (Peripheral Blood Mononuclear Cells) single-cell RNA-seq data and wants to annotate cell types.

## Input Data
```python
import scanpy as sc

# Load preprocessed PBMC data
adata = sc.read_h5ad("pbmc_5k.h5ad")

# Data characteristics:
# - 5,000 cells
# - 20,000 genes
# - Already normalized, log-transformed
# - Leiden clustering performed (15 clusters)
```

## Execution Walkthrough

### Phase 0: Data QC Check

**System checks**:
```
✓ Data format: h5ad (AnnData)
✓ Preprocessing: normalized, log1p applied
✓ Clustering: leiden (15 clusters)
✓ Gene names: standardized (uppercase)
✓ Quality metrics: median 1,200 genes/cell, 3% MT
```

**Decision**: PASS → Proceed to Phase 1

---

### Phase 1: Annotation Strategy Selection

**User input**:
- Tissue: PBMC
- Species: human
- No custom markers provided

**Agent 1: Tissue Biology Expert**

```yaml
Assessment:
  tissue_validation: PASS
  expected_cell_types:
    - name: "T cells"
      proportion: "60-70%"
      markers: [CD3D, CD3E, CD3G]
    - name: "B cells"
      proportion: "10-15%"
      markers: [CD19, MS4A1, CD79A]
    - name: "NK cells"
      proportion: "5-10%"
      markers: [NKG7, GNLY, NCAM1]
    - name: "Monocytes"
      proportion: "10-20%"
      markers: [CD14, LYZ, S100A8]
    - name: "Dendritic cells"
      proportion: "1-3%"
      markers: [FCER1A, CD1C, CLEC10A]

  clustering_assessment: OPTIMAL
  verdict: PASS
```

**Agent 2: Method Specialist**

```yaml
Recommendation:
  primary_method: "hybrid"
  strategy:
    step1:
      method: "marker-based"
      target: "major cell types"
    step2:
      method: "reference-based"
      tool: "SingleR"
      reference: "Monaco_Immune_2019"
      target: "subtype resolution"
    step3:
      method: "llm-assisted"
      target: "ambiguous clusters"

  confidence: HIGH
  verdict: PASS
```

**Agent 3: Quality Controller**

```yaml
Quality_Assessment:
  data_quality: GOOD
  clustering_quality: GOOD
  batch_effects: NONE
  doublet_risk: LOW (3 clusters flagged)

  recommendations:
    - "Run doublet detection on clusters 5, 11, 14"
    - "Cluster 8 has low gene count - flag for review"

  verdict: CONDITIONAL
  conditions:
    - "Address doublet concern for flagged clusters"
```

**Deliberation**: 2 rounds → Converged

**Output**: `annotation-strategy.yaml`

```yaml
primary_method: "hybrid"
expected_cell_types: [5 major types listed above]
strategy_details:
  step1:
    method: "marker-based"
    confidence_threshold: 0.7
  step2:
    method: "reference-based"
    tool: "SingleR"
    reference: "Monaco_Immune_2019"
  validation:
    cross_method_agreement_threshold: 0.9
```

**User Prompt**:
```
Phase 1 Complete: Annotation Strategy Approved

Strategy: Hybrid (marker-based + SingleR)
Expected cell types: 5 major types
Confidence: HIGH

⚠️  Note: 3 clusters flagged for doublet check

Proceed to Phase 2 (Marker Library Preparation)?
[Approve] [Modify] [Cancel]
```

**User**: Approve

---

### Phase 2: Marker Gene Library Preparation

**Step 1: Database Query** (parallel subagent)

```
Querying CellMarker... ✓ (retrieved 8 cell types, 45 markers)
Querying PanglaoDB... ✓ (retrieved 8 cell types, 52 markers)
Querying ImmGen... ✓ (retrieved 6 cell types, 38 markers)

Consolidated: 8 cell types, 87 unique markers
```

**Step 2: Detection Validation** (subagent)

```
Checking marker detection in dataset...

T cells: 8/10 markers detected (80%)
  Missing: CD3A, CD3B
  Resolved: CD20 → MS4A1 (synonym)

B cells: 6/8 markers detected (75%)
  Missing: PAX5, BANK1

NK cells: 5/6 markers detected (83%)
  Missing: CD56 (NCAM1 detected instead)

Monocytes: 7/7 markers detected (100%)

Dendritic cells: 3/5 markers detected (60%)
  Missing: CD1C, CLEC10A
  ⚠️  Low coverage - recommend reference-based for this type

Average detection rate: 78%
```

**Step 3: Prioritization & Conflict Resolution**

```
Ambiguous markers identified:
  - CD4: expressed in T cells AND monocytes
  - CD14: expressed in monocytes AND dendritic cells

Marker combinations created:
  - "CD4+ AND CD3+" → T cells
  - "CD4+ AND CD14+" → monocytes
  - "CD14+ AND HLA-DRA+" → monocytes
  - "CD14+ AND FCER1A+" → dendritic cells
```

**Output**: `marker_gene_library.yaml`

```yaml
cell_types:
  T cells:
    canonical_markers: [CD3D, CD3E, CD3G]
    auxiliary_markers: [CD2, CD5, CD7, CD8A, CD4]
    detection_rate: 0.80
    quality_grade: EXCELLENT

  B cells:
    canonical_markers: [CD19, MS4A1, CD79A]
    auxiliary_markers: [CD79B, BLK]
    detection_rate: 0.75
    quality_grade: GOOD

  # ... (other cell types)

  Dendritic cells:
    canonical_markers: [FCER1A, CLEC9A]
    auxiliary_markers: [CD1C, CLEC10A]
    detection_rate: 0.60
    quality_grade: MARGINAL
```

**User Prompt**:
```
Phase 2 Complete: Marker Library Ready

Summary:
- 8 cell types
- 87 markers (78% detected)
- 6 cell types: EXCELLENT/GOOD
- 1 cell type: MARGINAL (Dendritic cells)

⚠️  Dendritic cells: Only 60% markers detected
    → Will rely on reference-based annotation

Proceed to Phase 3 (Annotation Execution)?
[Approve] [Review Details] [Modify]
```

**User**: Approve

---

### Phase 3: Multi-Strategy Annotation Execution

**Parallel execution** (3 methods):

#### Method A: Marker-Based

```
Running marker-based annotation...

Cluster 0: T cells (score: 2.34, confidence: 0.85)
  Evidence: CD3D(3.2), CD3E(2.9), CD3G(2.1), CD2(1.8)

Cluster 1: T cells (score: 2.18, confidence: 0.82)
  Evidence: CD3D(2.8), CD3E(2.6), CD8A(3.1)

Cluster 2: B cells (score: 1.89, confidence: 0.72)
  Evidence: CD19(2.5), MS4A1(2.1), CD79A(1.8)

Cluster 3: NK cells (score: 2.45, confidence: 0.88)
  Evidence: NKG7(3.5), GNLY(3.2), NCAM1(2.1)

Cluster 4: Monocytes (score: 2.67, confidence: 0.91)
  Evidence: CD14(3.8), LYZ(3.5), S100A8(2.9)

Cluster 5: Dendritic cells (score: 1.12, confidence: 0.45)
  ⚠️  Low confidence - flagged for review

# ... (clusters 6-14)

Summary: 15/15 clusters annotated
High confidence: 10 clusters
Low confidence: 2 clusters (5, 11)
```

#### Method B: Reference-Based (SingleR)

```
Running SingleR with Monaco_Immune_2019...

Cluster 0: CD8+ T cells (confidence: 0.92)
Cluster 1: CD4+ T cells (confidence: 0.88)
Cluster 2: B cells (confidence: 0.85)
Cluster 3: NK cells (confidence: 0.90)
Cluster 4: Classical monocytes (confidence: 0.93)
Cluster 5: Dendritic cells (confidence: 0.68)

# ... (clusters 6-14)

Summary: 15/15 clusters annotated
Average confidence: 0.84
```

#### Method C: LLM-Assisted (for ambiguous clusters)

```
Running LLM annotation for clusters 5, 11...

Cluster 5:
  Top markers: FCER1A, HLA-DRA, CD1C, CLEC10A, ...

  LLM Response:
    Cell Type: Dendritic cells
    Confidence: MEDIUM
    Reasoning: Expression of FCER1A and HLA-DRA strongly suggests
    myeloid dendritic cells. CD1C is a canonical marker for cDC2
    dendritic cells. However, some monocyte markers are also present,
    suggesting this may be a transitional population.

Cluster 11:
  Top markers: CD3D, CD3E, CD4, IL7R, CCR7, ...

  LLM Response:
    Cell Type: Naive CD4+ T cells
    Confidence: HIGH
    Reasoning: Strong expression of CD3 complex (T cell lineage),
    CD4 (helper T cells), IL7R and CCR7 (naive T cell markers).
    Absence of memory markers (CD44, CD69) confirms naive state.
```

**Consensus Integration**:

```
Integrating results from 3 methods...

Cluster 0:
  Marker-based: T cells (0.85)
  SingleR: CD8+ T cells (0.92)
  Consensus: CD8+ T cells (0.89, STRONG_MAJORITY)

Cluster 5:
  Marker-based: Dendritic cells (0.45)
  SingleR: Dendritic cells (0.68)
  LLM: Dendritic cells (0.70)
  Consensus: Dendritic cells (0.61, WEAK_MAJORITY)
  ⚠️  Flagged for manual review

# ... (other clusters)

Summary:
- Full agreement: 10 clusters (66.7%)
- Strong majority: 3 clusters (20%)
- Weak majority: 2 clusters (13.3%)
- Flagged for review: 2 clusters
```

**Domain Sanity Check**:

```
Checking biological plausibility...

✓ T cells: 62% (expected 60-70%)
✓ B cells: 12% (expected 10-15%)
✓ NK cells: 8% (expected 5-10%)
✓ Monocytes: 15% (expected 10-20%)
⚠️  Dendritic cells: 2.5% (expected 1-3%) - slightly high but acceptable

✓ All proportions within expected ranges
✓ No impossible cell types detected
✓ Marker expression validates annotations

Overall: PASS
```

**User Prompt**:
```
Phase 3 Complete: Consensus Annotations Ready

Summary:
- 15 clusters annotated
- 10 full agreement (66.7%)
- Average confidence: 0.81
- 2 clusters flagged for review

Flagged Clusters:
⚠️  Cluster 5: Dendritic cells (confidence: 0.61)
    Methods: marker(0.45), SingleR(0.68), LLM(0.70)
    Recommendation: Accept with caution

⚠️  Cluster 11: Naive CD4+ T cells (confidence: 0.75)
    Methods: marker(0.72), SingleR(0.78), LLM(0.90)
    Recommendation: Accept

Proceed to Phase 4 (Results Integration)?
[Approve All] [Review Flagged] [Modify]
```

**User**: Approve All

---

### Phase 4: Results Integration & Validation

**Apply annotations**:

```
Applying consensus annotations to data...

✓ Added 'cell_type' column (5,000 cells)
✓ Added 'annotation_confidence' column
✓ Added 'agreement_status' column

Cell types identified:
  - CD8+ T cells: 1,850 cells (37%)
  - CD4+ T cells: 1,250 cells (25%)
  - B cells: 600 cells (12%)
  - NK cells: 400 cells (8%)
  - Classical monocytes: 550 cells (11%)
  - Non-classical monocytes: 200 cells (4%)
  - Dendritic cells: 125 cells (2.5%)
  - Plasma cells: 25 cells (0.5%)
```

**Generate figures**:

```
Generating publication-quality figures...

✓ umap_cell_types.pdf (UMAP colored by cell type)
✓ umap_confidence.pdf (UMAP colored by confidence)
✓ cell_type_proportions.pdf (bar plot)
✓ marker_heatmap.pdf (canonical markers × cell types)
✓ confidence_distribution.pdf (histogram)

All figures saved to docs/05_figures/
```

**Evidence report**:

```markdown
# Cell Type Annotation Evidence Report

## CD8+ T cells (1,850 cells, 37%)
- Average confidence: 0.89
- Marker validation: PASS
  - CD3D: 85% cells expressing (mean: 2.8)
  - CD3E: 82% cells expressing (mean: 2.6)
  - CD8A: 88% cells expressing (mean: 3.1)
- Method agreement: FULL (marker + SingleR)

## Dendritic cells (125 cells, 2.5%)
- Average confidence: 0.61
- Marker validation: PASS
  - FCER1A: 65% cells expressing (mean: 2.1)
  - HLA-DRA: 78% cells expressing (mean: 2.5)
- Method agreement: WEAK_MAJORITY
- Note: Low marker coverage in library, relied on reference-based
```

**User Prompt**:
```
Phase 4 Complete: Results Integrated

Summary:
- 5,000 cells annotated
- 8 cell types identified
- Average confidence: 0.81
- Validation: PASS

All figures and reports generated.

Proceed to Phase 5 (Final Export)?
[Approve] [Review Evidence] [Regenerate Figures]
```

**User**: Approve

---

### Phase 5: Verification & Export

**Reproducibility package**:

```
Generating reproducible pipeline...

✓ code/annotation_pipeline.py (standalone script)
✓ docs/06_results/environment.yaml
✓ docs/06_results/requirements.txt
✓ docs/06_results/README.md
```

**Quality metrics**:

```
Computing quality metrics...

Annotation Coverage:
  - Total cells: 5,000
  - Annotated: 5,000 (100%)
  - High confidence (≥0.8): 4,200 (84%)

Confidence Metrics:
  - Mean: 0.81
  - Median: 0.85
  - Min: 0.45
  - Max: 0.93

Marker Consistency:
  - CD8+ T cells: 85% avg marker expression
  - B cells: 78% avg marker expression
  - (all cell types validated)

Quality Grade: EXCELLENT
```

**Documentation**:

```
Generating publication materials...

✓ methods_section.txt (ready for paper)
✓ figure_legends.txt (for all 5 figures)
✓ table_s1_cell_type_summary.csv
✓ table_s2_marker_genes.csv
```

**Final package**:

```
Packaging deliverables...

✓ cell_annotation_results.tar.gz (12.3 MB)

Contents:
  - figures/ (5 PDFs)
  - results/ (annotated data + tables)
  - code/ (reproducible pipeline)
  - docs/ (methods + legends)
```

**Final Summary**:

```
═══════════════════════════════════════════════════════════
Cell Type Annotation Complete

Summary Statistics:
- Total cells: 5,000
- Cell types identified: 8
- Average confidence: 0.81
- High confidence cells (≥0.8): 4,200 (84%)

Quality Grade: EXCELLENT

Cell Type Distribution:
- CD8+ T cells: 1,850 cells (37%)
- CD4+ T cells: 1,250 cells (25%)
- B cells: 600 cells (12%)
- Classical monocytes: 550 cells (11%)
- NK cells: 400 cells (8%)
- Non-classical monocytes: 200 cells (4%)
- Dendritic cells: 125 cells (2.5%)
- Plasma cells: 25 cells (0.5%)

Validation Status: PASS

Deliverables:
All results saved to docs/05_figures/ and docs/06_results/
Complete package: cell_annotation_results.tar.gz

Next Steps:
1. Review evidence report: docs/06_results/evidence_report.md
2. Use annotated data: docs/06_results/annotated_data.h5ad
3. Cite methods: docs/06_results/methods_section.txt

Thank you for using Amplify Cell Annotation!
═══════════════════════════════════════════════════════════
```

---

## Total Execution Time

- Phase 0: 30 seconds
- Phase 1: 5 minutes (2 deliberation rounds)
- Phase 2: 3 minutes (database queries + validation)
- Phase 3: 8 minutes (3 methods in parallel)
- Phase 4: 2 minutes (figures + validation)
- Phase 5: 2 minutes (exports + documentation)

**Total: ~20 minutes** (with user approval pauses)

---

## Key Takeaways

1. **Multi-agent deliberation** caught potential doublet issues early
2. **Hybrid strategy** compensated for low marker coverage in dendritic cells
3. **Domain sanity check** validated proportions matched PBMC expectations
4. **Consensus voting** resolved ambiguities with high confidence
5. **Complete documentation** ready for publication

This example demonstrates the full power of the cell annotation skill: rigorous, evidence-based, biologically validated, and publication-ready.
