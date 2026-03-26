#!/usr/bin/env python3
"""
Bio-Aging-Clocks: Complete Example Pipeline
============================================
Demonstrates three workflows:
  1. DNAm clocks on GEO data
  2. Transcriptomic clocks on local RNA-seq data (e.g., TCGA)
  3. Multi-clock comparison with visualization

Usage:
  python run_aging_clocks.py --workflow dnam     # GEO methylation data
  python run_aging_clocks.py --workflow rna      # Local RNA-seq data
  python run_aging_clocks.py --workflow compare  # Multi-clock comparison
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd


# ============================================================
# Workflow 1: DNAm Clocks on GEO Data
# ============================================================
def workflow_dnam(geo_id="GSE30870"):
    """Run epigenetic aging clocks on a GEO methylation dataset."""
    from biolearn.data_library import DataLibrary
    from biolearn.model_gallery import ModelGallery

    print(f"[1/4] Loading GEO dataset {geo_id}...")
    data = DataLibrary().get(geo_id).load()
    print(f"       Matrix: {data.dnam.shape[0]} CpGs × {data.dnam.shape[1]} samples")
    print(f"       Metadata columns: {list(data.metadata.columns)}")

    # Quality check
    print("[2/4] Running quality assessment...")
    report = data.quality_report()
    print(f"       Sample report shape: {report.sample_report.shape}")

    # Run clock panel
    print("[3/4] Running aging clocks...")
    gallery = ModelGallery()
    clock_names = ["Horvathv1", "Hannum", "PhenoAge", "DNAmTL"]
    results = data.metadata[["age"]].copy()

    for name in clock_names:
        model = gallery.get(name)
        pred = model.predict(data)
        results[name] = pred["Predicted"]

        valid = results[["age", name]].dropna()
        r = valid["age"].corr(valid[name])
        mae = np.abs(valid["age"] - valid[name]).mean()
        print(f"       {name:20s}  r={r:.4f}  MAE={mae:.2f}")

    # Age acceleration
    print("[4/4] Computing age acceleration...")
    from sklearn.linear_model import LinearRegression

    for name in clock_names:
        valid = results[["age", name]].dropna()
        lr = LinearRegression().fit(valid[["age"]], valid[name])
        results[f"{name}_accel"] = results[name] - lr.predict(results[["age"]])

    output_path = f"aging_clocks_{geo_id}.csv"
    results.to_csv(output_path)
    print(f"       Results saved to {output_path}")
    return results


# ============================================================
# Workflow 2: Transcriptomic Clocks on Local RNA-seq Data
# ============================================================
def workflow_rna(counts_path, clinical_path, sep="\t"):
    """Run transcriptomic aging clocks on local RNA-seq data (e.g., TCGA).

    Args:
        counts_path: Path to expression matrix (genes × samples, Ensembl IDs)
        clinical_path: Path to clinical metadata (must have sample, age columns)
        sep: Delimiter for input files
    """
    import biolearn
    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery

    print(f"[1/4] Loading expression data from {counts_path}...")
    compression = "gzip" if counts_path.endswith(".gz") else None
    counts = pd.read_csv(counts_path, sep=sep, compression=compression)

    # Identify gene ID column (first non-numeric column)
    gene_col = counts.columns[0]
    counts[gene_col] = counts[gene_col].str.split(".").str[0]  # strip version
    counts = counts.set_index(gene_col).groupby(level=0).mean()
    print(f"       {counts.shape[0]} genes × {counts.shape[1]} samples")

    print(f"[2/4] Loading clinical data from {clinical_path}...")
    clin = pd.read_csv(clinical_path, sep=sep, compression=compression)
    # Auto-detect age column
    age_col = [c for c in clin.columns if "age" in c.lower()][0]
    sample_col = clin.columns[0]
    clin_sub = clin[[sample_col, age_col]].copy()
    clin_sub.columns = ["sample", "age"]
    clin_sub = clin_sub.dropna(subset=["age"]).set_index("sample")

    common = counts.columns.intersection(clin_sub.index)
    print(f"       {len(common)} samples with both expression and age data")

    # CRITICAL: use rna= parameter
    geo_data = GeoData(metadata=clin_sub.loc[common], rna=counts[common])

    # Workaround: chdir to site-packages for Pasta relative path bug
    print("[3/4] Running transcriptomic clocks (with workarounds)...")
    biolearn_parent = os.path.dirname(os.path.dirname(biolearn.__file__))
    original_cwd = os.getcwd()
    os.chdir(biolearn_parent)

    gallery = ModelGallery()
    results = clin_sub.loc[common, ["age"]].copy()
    clock_names = ["Pasta", "REG", "TranscriptomicPredictionModel"]

    for name in clock_names:
        try:
            # Workaround: bypass ImputationDecorator
            pred = gallery.get(name).clock.predict(geo_data)
            results[name] = pred["Predicted"]
            valid = results[["age", name]].dropna()
            r = valid["age"].corr(valid[name])
            mae = np.abs(valid["age"] - valid[name]).mean()
            print(f"       {name:35s}  r={r:.4f}  MAE={mae:.2f}")
        except Exception as e:
            print(f"       {name:35s}  FAILED: {e}")

    os.chdir(original_cwd)

    print("[4/4] Saving results...")
    output_path = "aging_clocks_transcriptomic.csv"
    results.to_csv(output_path)
    print(f"       Results saved to {output_path}")
    return results


# ============================================================
# Workflow 3: Multi-Clock Comparison with Visualization
# ============================================================
def workflow_compare(geo_id="GSE30870"):
    """Run multiple clocks and generate comparison visualizations."""
    import matplotlib.pyplot as plt
    from biolearn.data_library import DataLibrary
    from biolearn.model_gallery import ModelGallery

    print(f"[1/3] Loading {geo_id}...")
    data = DataLibrary().get(geo_id).load()
    gallery = ModelGallery()

    clock_names = ["Horvathv1", "Hannum", "PhenoAge", "DNAmTL"]
    results = data.metadata[["age"]].copy()

    print("[2/3] Running clocks...")
    for name in clock_names:
        pred = gallery.get(name).predict(data)
        results[name] = pred["Predicted"]

    print("[3/3] Generating visualizations...")
    n_clocks = len(clock_names)
    fig, axes = plt.subplots(1, n_clocks, figsize=(5 * n_clocks, 5))

    for ax, name in zip(axes, clock_names):
        valid = results[["age", name]].dropna()
        r = valid["age"].corr(valid[name])
        mae = np.abs(valid["age"] - valid[name]).mean()

        ax.scatter(valid["age"], valid[name], alpha=0.4, s=20, c="steelblue")
        lims = [
            min(valid["age"].min(), valid[name].min()),
            max(valid["age"].max(), valid[name].max()),
        ]
        ax.plot(lims, lims, "r--", alpha=0.5)
        ax.set_xlabel("Chronological Age")
        ax.set_ylabel("Predicted Age")
        ax.set_title(f"{name}\nr={r:.3f}, MAE={mae:.1f}")

    plt.suptitle(f"Aging Clocks Comparison — {geo_id}", fontweight="bold")
    plt.tight_layout()

    output_fig = f"aging_clocks_comparison_{geo_id}.png"
    plt.savefig(output_fig, dpi=150, bbox_inches="tight")
    print(f"       Figure saved to {output_fig}")

    # Correlation matrix
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    corr_matrix = results[clock_names].corr()
    im = ax2.imshow(corr_matrix, cmap="RdYlBu_r", vmin=0, vmax=1)
    ax2.set_xticks(range(n_clocks))
    ax2.set_yticks(range(n_clocks))
    ax2.set_xticklabels(clock_names, rotation=45, ha="right")
    ax2.set_yticklabels(clock_names)
    for i in range(n_clocks):
        for j in range(n_clocks):
            ax2.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}", ha="center", va="center")
    plt.colorbar(im, ax=ax2)
    ax2.set_title("Inter-Clock Correlation")
    plt.tight_layout()

    corr_fig = f"aging_clocks_correlation_{geo_id}.png"
    plt.savefig(corr_fig, dpi=150, bbox_inches="tight")
    print(f"       Correlation matrix saved to {corr_fig}")

    return results


# ============================================================
# CLI Entry Point
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bio-Aging-Clocks Pipeline")
    parser.add_argument(
        "--workflow",
        choices=["dnam", "rna", "compare"],
        default="dnam",
        help="Which workflow to run",
    )
    parser.add_argument("--geo-id", default="GSE30870", help="GEO accession ID")
    parser.add_argument("--counts", help="Path to expression matrix (for rna workflow)")
    parser.add_argument("--clinical", help="Path to clinical metadata (for rna workflow)")

    args = parser.parse_args()

    if args.workflow == "dnam":
        workflow_dnam(args.geo_id)
    elif args.workflow == "rna":
        if not args.counts or not args.clinical:
            print("ERROR: --counts and --clinical are required for rna workflow")
            sys.exit(1)
        workflow_rna(args.counts, args.clinical)
    elif args.workflow == "compare":
        workflow_compare(args.geo_id)
