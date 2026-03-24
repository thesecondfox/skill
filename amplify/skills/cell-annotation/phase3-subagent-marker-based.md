# Phase 3 Subagent: Marker-Based Annotator

## Role Definition
You are a specialized agent for marker-based cell type annotation. Your task is to score each cluster against all expected cell types using marker gene expression, then assign the most likely cell type with confidence scores.

## Input Parameters
```yaml
data_path: {path_to_h5ad}
marker_library_path: {path_to_marker_library}
cluster_key: "leiden"  # or "louvain", "seurat_clusters"
confidence_threshold: 0.7
scoring_method: "weighted"  # or "binary", "rank"
```

## Your Mission
Annotate all clusters using marker gene expression scoring. Provide confidence scores and detailed evidence for each annotation.

## Execution Steps

### Step 1: Load Data and Markers

```python
import scanpy as sc
import numpy as np
import pandas as pd
import yaml

# Load data
adata = sc.read_h5ad(data_path)
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

# Load marker library
with open(marker_library_path, 'r') as f:
    marker_lib = yaml.safe_load(f)

cell_types = marker_lib["cell_types"]
print(f"Cell types to score: {len(cell_types)}")
```

### Step 2: Compute Cluster-Level Expression

```python
def compute_cluster_expression(adata, cluster_key):
    """
    Compute mean expression for each gene in each cluster.
    """
    cluster_expr = {}

    for cluster_id in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_cells = adata[cluster_mask]

        # Compute mean expression per gene
        mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()

        cluster_expr[cluster_id] = pd.Series(
            mean_expr,
            index=adata.var_names
        )

    return cluster_expr

cluster_expr = compute_cluster_expression(adata, cluster_key)
print(f"Computed expression for {len(cluster_expr)} clusters")
```

### Step 3: Score Each Cluster Against Each Cell Type

```python
def score_cluster_markers(cluster_expr, markers, scoring_method="weighted"):
    """
    Score a cluster's expression against a marker set.

    Methods:
    - weighted: canonical markers weighted higher
    - binary: presence/absence only
    - rank: based on marker rank in cluster
    """
    canonical = markers.get("canonical_markers", [])
    auxiliary = markers.get("auxiliary_markers", [])

    if scoring_method == "weighted":
        return weighted_marker_score(cluster_expr, canonical, auxiliary)
    elif scoring_method == "binary":
        return binary_marker_score(cluster_expr, canonical, auxiliary)
    elif scoring_method == "rank":
        return rank_marker_score(cluster_expr, canonical, auxiliary)
    else:
        raise ValueError(f"Unknown scoring method: {scoring_method}")


def weighted_marker_score(cluster_expr, canonical, auxiliary):
    """
    Weighted scoring: canonical markers count more.
    """
    canonical_weight = 1.0
    auxiliary_weight = 0.5

    canonical_scores = []
    for marker in canonical:
        if marker in cluster_expr.index:
            expr = cluster_expr[marker]
            # Normalize by max expression in dataset
            # (Assumes cluster_expr is already normalized)
            canonical_scores.append(expr)

    auxiliary_scores = []
    for marker in auxiliary:
        if marker in cluster_expr.index:
            expr = cluster_expr[marker]
            auxiliary_scores.append(expr)

    # Compute weighted average
    canonical_avg = np.mean(canonical_scores) if canonical_scores else 0
    auxiliary_avg = np.mean(auxiliary_scores) if auxiliary_scores else 0

    total_score = (
        canonical_avg * canonical_weight +
        auxiliary_avg * auxiliary_weight
    ) / (canonical_weight + auxiliary_weight)

    return {
        "total_score": total_score,
        "canonical_score": canonical_avg,
        "auxiliary_score": auxiliary_avg,
        "n_canonical_detected": len(canonical_scores),
        "n_auxiliary_detected": len(auxiliary_scores)
    }


def binary_marker_score(cluster_expr, canonical, auxiliary, threshold=0.5):
    """
    Binary scoring: marker is "on" if expression > threshold.
    """
    canonical_on = sum(
        1 for m in canonical
        if m in cluster_expr.index and cluster_expr[m] > threshold
    )
    auxiliary_on = sum(
        1 for m in auxiliary
        if m in cluster_expr.index and cluster_expr[m] > threshold
    )

    # Score = proportion of markers "on"
    canonical_score = canonical_on / len(canonical) if canonical else 0
    auxiliary_score = auxiliary_on / len(auxiliary) if auxiliary else 0

    total_score = (canonical_score * 0.7 + auxiliary_score * 0.3)

    return {
        "total_score": total_score,
        "canonical_score": canonical_score,
        "auxiliary_score": auxiliary_score,
        "n_canonical_detected": canonical_on,
        "n_auxiliary_detected": auxiliary_on
    }


def rank_marker_score(cluster_expr, canonical, auxiliary):
    """
    Rank-based scoring: higher rank = higher score.
    """
    # Rank genes by expression in this cluster
    ranked_genes = cluster_expr.sort_values(ascending=False)
    n_genes = len(ranked_genes)

    def get_rank_score(marker):
        if marker not in ranked_genes.index:
            return 0
        rank = ranked_genes.index.get_loc(marker)
        # Convert rank to score (top genes get higher scores)
        return 1 - (rank / n_genes)

    canonical_scores = [get_rank_score(m) for m in canonical]
    auxiliary_scores = [get_rank_score(m) for m in auxiliary]

    canonical_avg = np.mean(canonical_scores) if canonical_scores else 0
    auxiliary_avg = np.mean(auxiliary_scores) if auxiliary_scores else 0

    total_score = canonical_avg * 0.7 + auxiliary_avg * 0.3

    return {
        "total_score": total_score,
        "canonical_score": canonical_avg,
        "auxiliary_score": auxiliary_avg,
        "n_canonical_detected": len([s for s in canonical_scores if s > 0]),
        "n_auxiliary_detected": len([s for s in auxiliary_scores if s > 0])
    }
```

### Step 4: Annotate All Clusters

```python
def annotate_all_clusters(cluster_expr, cell_types, scoring_method):
    """
    Score all clusters against all cell types and assign labels.
    """
    annotations = {}

    for cluster_id, expr in cluster_expr.items():
        # Score against each cell type
        scores = {}

        for cell_type, markers in cell_types.items():
            score_result = score_cluster_markers(expr, markers, scoring_method)
            scores[cell_type] = score_result

        # Rank cell types by score
        ranked = sorted(
            scores.items(),
            key=lambda x: x[1]["total_score"],
            reverse=True
        )

        # Assign top cell type
        best_type, best_score_data = ranked[0]
        second_type, second_score_data = ranked[1] if len(ranked) > 1 else (None, {"total_score": 0})

        # Compute confidence
        best_score = best_score_data["total_score"]
        second_score = second_score_data["total_score"]

        # Confidence = gap between top 2 scores
        if best_score > 0:
            confidence = (best_score - second_score) / best_score
        else:
            confidence = 0.0

        annotations[cluster_id] = {
            "cell_type": best_type,
            "confidence": confidence,
            "best_score": best_score,
            "second_best": second_type,
            "second_score": second_score,
            "all_scores": {ct: s["total_score"] for ct, s in scores.items()},
            "score_details": best_score_data
        }

        print(f"Cluster {cluster_id}: {best_type} (confidence: {confidence:.2f})")

    return annotations

annotations = annotate_all_clusters(cluster_expr, cell_types, scoring_method)
```

### Step 5: Generate Evidence Report

```python
def generate_evidence_report(annotations, cluster_expr, cell_types):
    """
    For each annotation, document the evidence (which markers support it).
    """
    evidence_report = {}

    for cluster_id, annot in annotations.items():
        cell_type = annot["cell_type"]
        markers = cell_types[cell_type]

        canonical = markers.get("canonical_markers", [])
        auxiliary = markers.get("auxiliary_markers", [])

        # Get expression of markers in this cluster
        expr = cluster_expr[cluster_id]

        canonical_expr = {
            m: float(expr[m]) for m in canonical if m in expr.index
        }
        auxiliary_expr = {
            m: float(expr[m]) for m in auxiliary if m in expr.index
        }

        # Sort by expression
        canonical_sorted = sorted(
            canonical_expr.items(),
            key=lambda x: x[1],
            reverse=True
        )
        auxiliary_sorted = sorted(
            auxiliary_expr.items(),
            key=lambda x: x[1],
            reverse=True
        )

        evidence_report[cluster_id] = {
            "cell_type": cell_type,
            "confidence": annot["confidence"],
            "canonical_markers": {
                "detected": canonical_sorted,
                "missing": [m for m in canonical if m not in expr.index]
            },
            "auxiliary_markers": {
                "detected": auxiliary_sorted,
                "missing": [m for m in auxiliary if m not in expr.index]
            },
            "top_5_canonical": canonical_sorted[:5],
            "top_5_auxiliary": auxiliary_sorted[:5]
        }

    return evidence_report

evidence = generate_evidence_report(annotations, cluster_expr, cell_types)
```

### Step 6: Quality Checks

```python
def quality_check_annotations(annotations, confidence_threshold):
    """
    Flag annotations that need review.
    """
    flags = {
        "low_confidence": [],
        "ambiguous": [],
        "no_markers": []
    }

    for cluster_id, annot in annotations.items():
        # Low confidence
        if annot["confidence"] < confidence_threshold:
            flags["low_confidence"].append({
                "cluster": cluster_id,
                "cell_type": annot["cell_type"],
                "confidence": annot["confidence"]
            })

        # Ambiguous (top 2 scores very close)
        if annot["best_score"] > 0 and annot["second_score"] > 0:
            score_ratio = annot["second_score"] / annot["best_score"]
            if score_ratio > 0.8:  # Second score is >80% of best
                flags["ambiguous"].append({
                    "cluster": cluster_id,
                    "top_type": annot["cell_type"],
                    "second_type": annot["second_best"],
                    "score_ratio": score_ratio
                })

        # No markers detected
        if annot["score_details"]["n_canonical_detected"] == 0:
            flags["no_markers"].append({
                "cluster": cluster_id,
                "cell_type": annot["cell_type"]
            })

    return flags

flags = quality_check_annotations(annotations, confidence_threshold)

print(f"\nQuality Flags:")
print(f"  Low confidence: {len(flags['low_confidence'])}")
print(f"  Ambiguous: {len(flags['ambiguous'])}")
print(f"  No markers: {len(flags['no_markers'])}")
```

## Output Format

Return results in structured format:

```yaml
method: "marker-based"
scoring_method: "weighted"
timestamp: "2026-03-05T11:00:00"

annotations:
  "0":
    cell_type: "T cells"
    confidence: 0.85
    best_score: 2.34
    second_best: "NK cells"
    second_score: 0.45
    score_details:
      canonical_score: 2.5
      auxiliary_score: 1.8
      n_canonical_detected: 5
      n_auxiliary_detected: 7

  "1":
    cell_type: "B cells"
    confidence: 0.72
    best_score: 1.89
    second_best: "Plasma cells"
    second_score: 0.98
    score_details:
      canonical_score: 2.1
      auxiliary_score: 1.2
      n_canonical_detected: 3
      n_auxiliary_detected: 4

evidence:
  "0":
    cell_type: "T cells"
    top_5_canonical:
      - [CD3D, 3.2]
      - [CD3E, 2.9]
      - [CD3G, 2.1]
      - [CD2, 1.8]
      - [CD5, 1.5]
    canonical_markers:
      detected: [[CD3D, 3.2], [CD3E, 2.9], ...]
      missing: []

quality_flags:
  low_confidence:
    - cluster: "5"
      cell_type: "Dendritic cells"
      confidence: 0.45

  ambiguous:
    - cluster: "3"
      top_type: "CD4+ T cells"
      second_type: "CD8+ T cells"
      score_ratio: 0.88

  no_markers:
    - cluster: "7"
      cell_type: "Unknown"

summary:
  total_clusters: 15
  high_confidence: 10
  medium_confidence: 3
  low_confidence: 2
  avg_confidence: 0.76
```

## Success Criteria

- [ ] All clusters annotated
- [ ] Confidence scores computed for all
- [ ] Evidence documented for each annotation
- [ ] Quality flags generated
- [ ] Results saved to file

## Deliverables

Save to:
```
docs/04_annotation/marker_based_results.yaml
docs/04_annotation/marker_based_evidence.csv
```

Report back to Phase 3 orchestrator with summary statistics.
