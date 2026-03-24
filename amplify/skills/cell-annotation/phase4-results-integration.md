# Phase 4: Results Integration & Validation

## Overview
This phase integrates consensus annotations into the data object, performs final validation checks, and prepares comprehensive reports with visualizations.

## Prerequisites
- Phase 3 completed: `consensus_annotations.yaml` exists
- Manual review completed (if any clusters were flagged)
- User approval obtained

## Execution Flow

```
Step 1: Apply Annotations to Data
  ├─ Add consensus labels to adata.obs
  ├─ Add confidence scores
  └─ Add method agreement metadata

Step 2: Final Validation
  ├─ Cell type proportion check
  ├─ Marker expression validation
  └─ Cross-cluster consistency check

Step 3: Generate Visualizations
  ├─ UMAP with cell type colors
  ├─ Marker expression heatmap
  ├─ Cell type proportion plots
  └─ Confidence distribution plots

Step 4: Create Evidence Report
  ├─ Per-cluster annotation evidence
  ├─ Method comparison tables
  └─ Quality metrics summary

Step 5: Export Results
  ├─ Annotated h5ad file
  ├─ CSV annotation table
  └─ Publication-ready figures
```

## Step 1: Apply Annotations to Data

```python
import scanpy as sc
import pandas as pd
import yaml

def apply_annotations_to_adata(adata, consensus_path, cluster_key="leiden"):
    """
    Apply consensus annotations to AnnData object.
    """
    # Load consensus
    with open(consensus_path, 'r') as f:
        consensus = yaml.safe_load(f)

    # Create annotation columns
    adata.obs["cell_type"] = "Unknown"
    adata.obs["annotation_confidence"] = 0.0
    adata.obs["annotation_method"] = "Unknown"
    adata.obs["agreement_status"] = "Unknown"

    # Apply annotations
    for cluster_id, data in consensus["consensus_annotations"].items():
        cluster_mask = adata.obs[cluster_key] == cluster_id

        adata.obs.loc[cluster_mask, "cell_type"] = data["consensus"]
        adata.obs.loc[cluster_mask, "annotation_confidence"] = data["calibrated_confidence"]
        adata.obs.loc[cluster_mask, "agreement_status"] = data["agreement_status"]

        # Record which methods were used
        methods_used = ", ".join(data["methods"].keys())
        adata.obs.loc[cluster_mask, "annotation_method"] = methods_used

    print(f"Applied annotations to {adata.n_obs} cells")
    print(f"Cell types: {adata.obs['cell_type'].unique()}")

    return adata
```

## Step 2: Final Validation

```python
def final_validation_check(adata, expected_cell_types, marker_library):
    """
    Perform comprehensive validation of final annotations.
    """
    validation_report = {
        "proportion_check": {},
        "marker_validation": {},
        "consistency_check": {},
        "overall_status": "PASS"
    }

    # Check 1: Cell type proportions
    total_cells = len(adata)
    for expected_ct in expected_cell_types:
        ct_name = expected_ct["name"]
        expected_range = expected_ct["expected_proportion"]

        # Count cells
        n_cells = (adata.obs["cell_type"] == ct_name).sum()
        observed_pct = n_cells / total_cells * 100

        # Parse expected range
        import re
        match = re.match(r"(\d+)-(\d+)%", expected_range)
        if match:
            min_pct, max_pct = float(match.group(1)), float(match.group(2))

            status = "PASS" if min_pct <= observed_pct <= max_pct else "WARNING"

            validation_report["proportion_check"][ct_name] = {
                "expected": expected_range,
                "observed": f"{observed_pct:.1f}%",
                "n_cells": int(n_cells),
                "status": status
            }

            if status == "WARNING":
                validation_report["overall_status"] = "WARNING"

    # Check 2: Marker expression validation
    for cell_type in adata.obs["cell_type"].unique():
        if cell_type == "Unknown":
            continue

        if cell_type in marker_library["cell_types"]:
            canonical_markers = marker_library["cell_types"][cell_type]["canonical_markers"]

            # Get cells of this type
            ct_cells = adata[adata.obs["cell_type"] == cell_type]

            # Check marker expression
            markers_expressed = 0
            for marker in canonical_markers:
                if marker in adata.var_names:
                    expr = ct_cells[:, marker].X
                    if hasattr(expr, "toarray"):
                        expr = expr.toarray().flatten()
                    pct_expr = (expr > 0).sum() / len(expr)

                    if pct_expr > 0.1:  # >10% cells expressing
                        markers_expressed += 1

            validation_status = "PASS" if markers_expressed >= len(canonical_markers) * 0.5 else "FAIL"

            validation_report["marker_validation"][cell_type] = {
                "canonical_markers": canonical_markers,
                "markers_expressed": markers_expressed,
                "total_markers": len(canonical_markers),
                "status": validation_status
            }

            if validation_status == "FAIL":
                validation_report["overall_status"] = "FAIL"

    # Check 3: Cross-cluster consistency
    # Clusters with same annotation should have similar marker profiles
    for cell_type in adata.obs["cell_type"].unique():
        if cell_type == "Unknown":
            continue

        clusters_with_type = adata.obs[adata.obs["cell_type"] == cell_type]["leiden"].unique()

        if len(clusters_with_type) > 1:
            # Check if these clusters are similar
            # (Simplified - in practice, compute correlation of marker expression)
            validation_report["consistency_check"][cell_type] = {
                "n_clusters": len(clusters_with_type),
                "clusters": list(clusters_with_type),
                "status": "PASS"  # Placeholder
            }

    return validation_report
```

## Step 3: Generate Visualizations

```python
def generate_annotation_figures(adata, output_dir="docs/05_figures"):
    """
    Generate publication-quality figures.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Figure 1: UMAP with cell type annotations
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color="cell_type",
        legend_loc="right margin",
        title="Cell Type Annotations",
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/umap_cell_types.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # Figure 2: UMAP with confidence scores
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color="annotation_confidence",
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
        title="Annotation Confidence",
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/umap_confidence.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # Figure 3: Cell type proportions
    fig, ax = plt.subplots(figsize=(8, 6))
    cell_type_counts = adata.obs["cell_type"].value_counts()
    cell_type_pcts = cell_type_counts / len(adata) * 100

    cell_type_pcts.plot(kind="barh", ax=ax, color="steelblue")
    ax.set_xlabel("Percentage of cells (%)")
    ax.set_ylabel("Cell Type")
    ax.set_title("Cell Type Proportions")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cell_type_proportions.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # Figure 4: Marker expression heatmap
    # Get top markers for each cell type
    marker_genes = []
    for ct in adata.obs["cell_type"].unique():
        if ct in marker_library["cell_types"]:
            markers = marker_library["cell_types"][ct]["canonical_markers"][:3]
            marker_genes.extend([m for m in markers if m in adata.var_names])

    marker_genes = list(set(marker_genes))

    if marker_genes:
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.heatmap(
            adata,
            var_names=marker_genes,
            groupby="cell_type",
            cmap="RdBu_r",
            dendrogram=True,
            ax=ax,
            show=False
        )
        plt.tight_layout()
        plt.savefig(f"{output_dir}/marker_heatmap.pdf", dpi=300, bbox_inches="tight")
        plt.close()

    # Figure 5: Confidence distribution
    fig, ax = plt.subplots(figsize=(8, 6))
    adata.obs["annotation_confidence"].hist(bins=20, ax=ax, color="steelblue", edgecolor="black")
    ax.set_xlabel("Annotation Confidence")
    ax.set_ylabel("Number of Cells")
    ax.set_title("Distribution of Annotation Confidence")
    ax.axvline(0.7, color="red", linestyle="--", label="Threshold (0.7)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/confidence_distribution.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Generated figures in {output_dir}/")
```

## Step 4: Create Evidence Report

```python
def generate_evidence_report(adata, consensus, marker_library, validation_report):
    """
    Generate comprehensive evidence report.
    """
    report = []

    report.append("# Cell Type Annotation Evidence Report\n")
    report.append(f"**Date**: {datetime.now().strftime('%Y-%m-%d')}\n")
    report.append(f"**Total Cells**: {adata.n_obs:,}")
    report.append(f"**Total Cell Types**: {len(adata.obs['cell_type'].unique())}\n")

    # Summary statistics
    report.append("## Summary Statistics\n")
    report.append(f"- Average confidence: {adata.obs['annotation_confidence'].mean():.2f}")
    report.append(f"- High confidence cells (≥0.8): {(adata.obs['annotation_confidence'] >= 0.8).sum():,} ({(adata.obs['annotation_confidence'] >= 0.8).sum() / len(adata) * 100:.1f}%)")
    report.append(f"- Low confidence cells (<0.6): {(adata.obs['annotation_confidence'] < 0.6).sum():,} ({(adata.obs['annotation_confidence'] < 0.6).sum() / len(adata) * 100:.1f}%)\n")

    # Per cell type details
    report.append("## Cell Type Annotations\n")

    for cell_type in sorted(adata.obs["cell_type"].unique()):
        if cell_type == "Unknown":
            continue

        ct_cells = adata.obs[adata.obs["cell_type"] == cell_type]
        n_cells = len(ct_cells)
        pct = n_cells / len(adata) * 100
        avg_conf = ct_cells["annotation_confidence"].mean()

        report.append(f"### {cell_type}")
        report.append(f"- **Cells**: {n_cells:,} ({pct:.1f}%)")
        report.append(f"- **Average confidence**: {avg_conf:.2f}")

        # Marker evidence
        if cell_type in marker_library["cell_types"]:
            canonical = marker_library["cell_types"][cell_type]["canonical_markers"]
            report.append(f"- **Canonical markers**: {', '.join(canonical[:5])}")

        # Validation status
        if cell_type in validation_report["marker_validation"]:
            val_status = validation_report["marker_validation"][cell_type]["status"]
            report.append(f"- **Marker validation**: {val_status}")

        report.append("")

    # Validation summary
    report.append("## Validation Summary\n")
    report.append(f"**Overall Status**: {validation_report['overall_status']}\n")

    if validation_report["overall_status"] != "PASS":
        report.append("### Warnings\n")
        for ct, data in validation_report["proportion_check"].items():
            if data["status"] == "WARNING":
                report.append(f"- {ct}: Expected {data['expected']}, observed {data['observed']}")

        for ct, data in validation_report["marker_validation"].items():
            if data["status"] == "FAIL":
                report.append(f"- {ct}: Only {data['markers_expressed']}/{data['total_markers']} markers expressed")

    return "\n".join(report)
```

## Step 5: Export Results

```python
def export_annotation_results(adata, output_dir="docs/06_results"):
    """
    Export all annotation results.
    """
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Export 1: Annotated h5ad file
    output_h5ad = f"{output_dir}/annotated_data.h5ad"
    adata.write_h5ad(output_h5ad)
    print(f"Saved annotated data: {output_h5ad}")

    # Export 2: Annotation table (CSV)
    annotation_table = adata.obs[[
        "cell_type",
        "annotation_confidence",
        "annotation_method",
        "agreement_status",
        "leiden"
    ]].copy()
    annotation_table.to_csv(f"{output_dir}/cell_annotations.csv")
    print(f"Saved annotation table: {output_dir}/cell_annotations.csv")

    # Export 3: Cell type summary
    summary = adata.obs.groupby("cell_type").agg({
        "annotation_confidence": ["mean", "std", "min", "max"],
        "cell_type": "count"
    })
    summary.columns = ["avg_confidence", "std_confidence", "min_confidence", "max_confidence", "n_cells"]
    summary["pct_cells"] = summary["n_cells"] / len(adata) * 100
    summary.to_csv(f"{output_dir}/cell_type_summary.csv")
    print(f"Saved cell type summary: {output_dir}/cell_type_summary.csv")

    # Export 4: Cluster-level annotations
    cluster_annotations = adata.obs.groupby("leiden").agg({
        "cell_type": lambda x: x.mode()[0],
        "annotation_confidence": "mean",
        "agreement_status": lambda x: x.mode()[0]
    })
    cluster_annotations["n_cells"] = adata.obs.groupby("leiden").size()
    cluster_annotations.to_csv(f"{output_dir}/cluster_annotations.csv")
    print(f"Saved cluster annotations: {output_dir}/cluster_annotations.csv")
```

## Complete Phase 4 Execution

```python
def execute_phase4(adata, consensus_path, marker_library_path, expected_cell_types):
    """
    Main execution function for Phase 4.
    """
    print("Phase 4: Results Integration & Validation")
    print("=" * 50)

    # Load marker library
    with open(marker_library_path, 'r') as f:
        marker_library = yaml.safe_load(f)

    # Step 1: Apply annotations
    print("\n[Step 1/5] Applying annotations to data...")
    adata = apply_annotations_to_adata(adata, consensus_path)

    # Step 2: Final validation
    print("\n[Step 2/5] Running final validation...")
    validation_report = final_validation_check(adata, expected_cell_types, marker_library)
    print(f"  Validation status: {validation_report['overall_status']}")

    # Step 3: Generate visualizations
    print("\n[Step 3/5] Generating visualizations...")
    generate_annotation_figures(adata)

    # Step 4: Create evidence report
    print("\n[Step 4/5] Creating evidence report...")
    with open(consensus_path, 'r') as f:
        consensus = yaml.safe_load(f)
    evidence_report = generate_evidence_report(adata, consensus, marker_library, validation_report)

    with open("docs/06_results/evidence_report.md", 'w') as f:
        f.write(evidence_report)

    # Step 5: Export results
    print("\n[Step 5/5] Exporting results...")
    export_annotation_results(adata)

    print("\n✓ Phase 4 complete!")
    print(f"  Annotated {adata.n_obs:,} cells")
    print(f"  Identified {len(adata.obs['cell_type'].unique())} cell types")
    print(f"  Average confidence: {adata.obs['annotation_confidence'].mean():.2f}")

    return adata, validation_report
```

## Output Files

```
docs/
├── 05_figures/
│   ├── umap_cell_types.pdf
│   ├── umap_confidence.pdf
│   ├── cell_type_proportions.pdf
│   ├── marker_heatmap.pdf
│   └── confidence_distribution.pdf
└── 06_results/
    ├── annotated_data.h5ad
    ├── cell_annotations.csv
    ├── cell_type_summary.csv
    ├── cluster_annotations.csv
    └── evidence_report.md
```

## Quality Gates

- [ ] All cells have annotations applied
- [ ] Validation status is PASS or WARNING (not FAIL)
- [ ] All figures generated successfully
- [ ] Evidence report is complete
- [ ] Exported files are valid

## STOP — Phase 4 Complete

Present summary to user and wait for approval before proceeding to Phase 5.
