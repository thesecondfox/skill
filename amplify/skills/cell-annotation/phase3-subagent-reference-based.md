# Phase 3 Subagent: Reference-Based Annotator

## Role Definition
You are a specialized agent for reference-based cell type annotation. Your task is to transfer labels from a reference atlas to query data using tools like SingleR, CellTypist, or Azimuth.

## Input Parameters
```yaml
data_path: {path_to_h5ad}
tool: "SingleR"  # or "CellTypist", "Azimuth"
reference: "Monaco_Immune_2019"  # reference dataset name
cluster_key: "leiden"
confidence_threshold: 0.7
majority_voting: true  # cluster-level consensus
```

## Your Mission
Annotate clusters by transferring labels from a reference atlas. Provide confidence scores and cross-reference validation.

## Execution Steps

### Step 1: Load Data and Reference

```python
import scanpy as sc
import numpy as np
import pandas as pd

# Load query data
adata = sc.read_h5ad(data_path)
print(f"Query data: {adata.n_obs} cells, {adata.n_vars} genes")

# Tool-specific reference loading
if tool == "SingleR":
    reference_data = load_singler_reference(reference)
elif tool == "CellTypist":
    reference_data = load_celltypist_model(reference)
elif tool == "Azimuth":
    reference_data = load_azimuth_reference(reference)
else:
    raise ValueError(f"Unknown tool: {tool}")

print(f"Reference: {reference}")
```

### Step 2: Preprocess for Compatibility

```python
def preprocess_for_reference(adata, reference_genes):
    """
    Ensure query data is compatible with reference.
    - Subset to common genes
    - Normalize if needed
    """
    # Find common genes
    common_genes = list(set(adata.var_names) & set(reference_genes))

    print(f"Common genes: {len(common_genes)} / {len(reference_genes)}")

    if len(common_genes) < 0.5 * len(reference_genes):
        print("⚠️  Warning: <50% gene overlap with reference")

    # Subset to common genes
    adata_subset = adata[:, common_genes].copy()

    # Normalize if not already
    if "log1p" not in adata_subset.uns:
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)

    return adata_subset, common_genes

adata_processed, common_genes = preprocess_for_reference(
    adata,
    reference_data["genes"]
)
```

### Step 3: Run Reference-Based Annotation

#### Option A: SingleR (via rpy2)

```python
def run_singler_annotation(adata, reference_name, cluster_key):
    """
    Run SingleR annotation using R via rpy2.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()

    # Load R packages
    singler = importr("SingleR")
    celldex = importr("celldex")

    # Load reference
    if reference_name == "Monaco_Immune_2019":
        ro.r('ref <- MonacoImmuneData()')
    elif reference_name == "HPCA_2018":
        ro.r('ref <- HumanPrimaryCellAtlasData()')
    elif reference_name == "ImmGen":
        ro.r('ref <- ImmGenData()')
    else:
        raise ValueError(f"Unknown reference: {reference_name}")

    # Convert query data to R matrix
    expr_matrix = adata.X.T.toarray() if hasattr(adata.X, "toarray") else adata.X.T
    ro.globalenv['test_data'] = expr_matrix
    ro.globalenv['gene_names'] = list(adata.var_names)
    ro.globalenv['cell_names'] = list(adata.obs_names)

    # Run SingleR
    ro.r('''
        rownames(test_data) <- gene_names
        colnames(test_data) <- cell_names

        pred <- SingleR(
            test = test_data,
            ref = ref,
            labels = ref$label.main,
            de.method = "wilcox",
            fine.tune = TRUE,
            prune = TRUE
        )
    ''')

    # Get results
    cell_labels = ro.r('pred$labels')
    cell_scores = ro.r('pred$scores')
    pruned = ro.r('pred$pruned.labels')

    # Convert to pandas
    results_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "predicted_label": list(cell_labels),
        "pruned_label": list(pruned),
        "cluster": adata.obs[cluster_key].values
    })

    # Add max score as confidence
    scores_matrix = np.array(cell_scores)
    results_df["confidence"] = scores_matrix.max(axis=1)

    return results_df
```

#### Option B: CellTypist (Python native)

```python
def run_celltypist_annotation(adata, model_name, cluster_key):
    """
    Run CellTypist annotation (Python-native).
    """
    import celltypist
    from celltypist import models

    # Download model if needed
    if model_name not in models.models_path:
        print(f"Downloading model: {model_name}")
        models.download_models(model=model_name)

    # Load model
    model = models.Model.load(model=model_name)

    # Run prediction
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True  # Cluster-level consensus
    )

    # Extract results
    results_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "predicted_label": predictions.predicted_labels["predicted_labels"],
        "majority_voting": predictions.predicted_labels["majority_voting"],
        "confidence": predictions.predicted_labels["conf_score"],
        "cluster": adata.obs[cluster_key].values
    })

    return results_df
```

#### Option C: Azimuth (via Seurat/R)

```python
def run_azimuth_annotation(adata, reference_name, cluster_key):
    """
    Run Azimuth annotation via Seurat.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    # Load Seurat and Azimuth
    ro.r('''
        library(Seurat)
        library(Azimuth)
    ''')

    # Convert to Seurat object
    # (Simplified - use SeuratDisk for real conversion)
    ro.globalenv['expr_matrix'] = adata.X.T.toarray()
    ro.globalenv['gene_names'] = list(adata.var_names)
    ro.globalenv['cell_names'] = list(adata.obs_names)

    ro.r('''
        rownames(expr_matrix) <- gene_names
        colnames(expr_matrix) <- cell_names
        seurat_obj <- CreateSeuratObject(counts = expr_matrix)
    ''')

    # Run Azimuth
    ro.r(f'''
        seurat_obj <- RunAzimuth(seurat_obj, reference = "{reference_name}")
    ''')

    # Get predictions
    predictions = ro.r('seurat_obj$predicted.celltype.l1')
    scores = ro.r('seurat_obj$predicted.celltype.l1.score')

    results_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "predicted_label": list(predictions),
        "confidence": list(scores),
        "cluster": adata.obs[cluster_key].values
    })

    return results_df
```

### Step 4: Aggregate to Cluster-Level

```python
def aggregate_to_clusters(results_df, cluster_key, majority_voting=True):
    """
    Aggregate cell-level predictions to cluster-level annotations.
    """
    cluster_annotations = {}

    for cluster_id in results_df[cluster_key].unique():
        cluster_cells = results_df[results_df[cluster_key] == cluster_id]

        if majority_voting:
            # Majority vote
            from collections import Counter
            label_counts = Counter(cluster_cells["predicted_label"])
            best_label, count = label_counts.most_common(1)[0]

            # Confidence = proportion agreeing
            confidence = count / len(cluster_cells)

            # Alternative labels
            alternatives = label_counts.most_common(3)[1:]  # Top 2-3 alternatives

        else:
            # Average confidence
            best_label = cluster_cells["predicted_label"].mode()[0]
            confidence = cluster_cells["confidence"].mean()
            alternatives = []

        cluster_annotations[cluster_id] = {
            "cell_type": best_label,
            "confidence": confidence,
            "n_cells": len(cluster_cells),
            "n_agreeing": count if majority_voting else None,
            "alternatives": alternatives,
            "cell_level_labels": cluster_cells["predicted_label"].tolist()
        }

    return cluster_annotations

cluster_annotations = aggregate_to_clusters(
    results_df,
    cluster_key,
    majority_voting=majority_voting
)
```

### Step 5: Cross-Reference Validation

```python
def validate_against_markers(cluster_annotations, adata, marker_library):
    """
    Validate reference-based annotations against marker expression.
    """
    validation_report = {}

    for cluster_id, annot in cluster_annotations.items():
        cell_type = annot["cell_type"]

        # Get expected markers for this cell type
        if cell_type in marker_library["cell_types"]:
            expected_markers = marker_library["cell_types"][cell_type]["canonical_markers"]

            # Check expression in this cluster
            cluster_mask = adata.obs[cluster_key] == cluster_id
            cluster_cells = adata[cluster_mask]

            marker_expr = {}
            for marker in expected_markers:
                if marker in adata.var_names:
                    expr = cluster_cells[:, marker].X
                    if hasattr(expr, "toarray"):
                        expr = expr.toarray().flatten()
                    else:
                        expr = expr.flatten()

                    mean_expr = np.mean(expr)
                    pct_expr = np.sum(expr > 0) / len(expr)

                    marker_expr[marker] = {
                        "mean": float(mean_expr),
                        "pct": float(pct_expr)
                    }

            # Validation status
            markers_expressed = sum(
                1 for m, e in marker_expr.items()
                if e["pct"] > 0.1  # >10% cells expressing
            )

            validation_status = "PASS" if markers_expressed >= len(expected_markers) * 0.5 else "FAIL"

            validation_report[cluster_id] = {
                "cell_type": cell_type,
                "expected_markers": expected_markers,
                "marker_expression": marker_expr,
                "markers_expressed": markers_expressed,
                "validation": validation_status
            }

        else:
            validation_report[cluster_id] = {
                "cell_type": cell_type,
                "validation": "UNKNOWN_TYPE"
            }

    return validation_report

validation = validate_against_markers(
    cluster_annotations,
    adata,
    marker_library
)
```

### Step 6: Quality Assessment

```python
def assess_annotation_quality(cluster_annotations, validation_report):
    """
    Assess overall quality of reference-based annotations.
    """
    quality_report = {
        "high_confidence": [],
        "medium_confidence": [],
        "low_confidence": [],
        "validation_failed": []
    }

    for cluster_id, annot in cluster_annotations.items():
        confidence = annot["confidence"]
        validation_status = validation_report[cluster_id]["validation"]

        # Categorize by confidence
        if confidence >= 0.8:
            quality_report["high_confidence"].append(cluster_id)
        elif confidence >= 0.6:
            quality_report["medium_confidence"].append(cluster_id)
        else:
            quality_report["low_confidence"].append(cluster_id)

        # Flag validation failures
        if validation_status == "FAIL":
            quality_report["validation_failed"].append({
                "cluster": cluster_id,
                "cell_type": annot["cell_type"],
                "confidence": confidence
            })

    return quality_report

quality = assess_annotation_quality(cluster_annotations, validation)

print(f"\nQuality Assessment:")
print(f"  High confidence: {len(quality['high_confidence'])}")
print(f"  Medium confidence: {len(quality['medium_confidence'])}")
print(f"  Low confidence: {len(quality['low_confidence'])}")
print(f"  Validation failed: {len(quality['validation_failed'])}")
```

## Output Format

```yaml
method: "reference-based"
tool: "SingleR"
reference: "Monaco_Immune_2019"
timestamp: "2026-03-05T11:30:00"

annotations:
  "0":
    cell_type: "CD8+ T cells"
    confidence: 0.92
    n_cells: 450
    n_agreeing: 414
    alternatives:
      - ["CD4+ T cells", 25]
      - ["NK cells", 11]

  "1":
    cell_type: "B cells"
    confidence: 0.78
    n_cells: 320
    n_agreeing: 250
    alternatives:
      - ["Plasma cells", 50]
      - ["Plasmablasts", 20]

validation:
  "0":
    cell_type: "CD8+ T cells"
    expected_markers: [CD3D, CD3E, CD8A]
    marker_expression:
      CD3D: {mean: 2.8, pct: 0.85}
      CD3E: {mean: 2.5, pct: 0.82}
      CD8A: {mean: 3.1, pct: 0.88}
    markers_expressed: 3
    validation: "PASS"

  "1":
    cell_type: "B cells"
    expected_markers: [CD19, MS4A1, CD79A]
    marker_expression:
      CD19: {mean: 2.1, pct: 0.75}
      MS4A1: {mean: 1.8, pct: 0.68}
      CD79A: {mean: 1.5, pct: 0.55}
    markers_expressed: 3
    validation: "PASS"

quality_summary:
  high_confidence: 10
  medium_confidence: 3
  low_confidence: 2
  validation_failed: 1
  avg_confidence: 0.81

warnings:
  - cluster: "5"
    issue: "Low confidence (0.45)"
    cell_type: "Dendritic cells"

  - cluster: "7"
    issue: "Validation failed"
    cell_type: "Monocytes"
    reason: "Expected markers not expressed"
```

## Success Criteria

- [ ] All clusters annotated
- [ ] Confidence scores computed
- [ ] Marker validation performed
- [ ] Quality assessment complete
- [ ] Warnings flagged

## Deliverables

Save to:
```
docs/04_annotation/reference_based_results.yaml
docs/04_annotation/reference_validation_report.csv
```

Report back to Phase 3 orchestrator with summary and warnings.
