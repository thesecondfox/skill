"""Preprocessing pipeline for perturbation benchmark data."""
import scanpy as sc
import numpy as np
from anndata import AnnData
from typing import List, Tuple, Optional


def preprocess_adata(
    adata: AnnData,
    condition_key: str = "condition",
    cell_type_key: str = "cell_type",
    control_label: str = "control",
    stim_label: str = "stimulated",
    n_top_genes: int = 2000,
    normalize: bool = True,
) -> AnnData:
    """Standard preprocessing pipeline for perturbation prediction.

    Steps:
      1. Filter genes/cells
      2. Normalize + log1p
      3. Select highly variable genes
      4. Validate required obs columns
    """
    adata = adata.copy()

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Normalize
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3",
                                 layer="counts")
    adata = adata[:, adata.var.highly_variable].copy()

    # Validate columns
    assert condition_key in adata.obs.columns, f"Missing column: {condition_key}"
    assert cell_type_key in adata.obs.columns, f"Missing column: {cell_type_key}"

    conditions = adata.obs[condition_key].unique()
    assert control_label in conditions, f"Control label '{control_label}' not found"
    assert stim_label in conditions, f"Stim label '{stim_label}' not found"

    return adata


def split_train_test(
    adata: AnnData,
    train_cell_types: List[str],
    test_cell_type: str,
    condition_key: str = "condition",
    cell_type_key: str = "cell_type",
    control_label: str = "control",
    stim_label: str = "stimulated",
) -> Tuple[AnnData, AnnData, AnnData]:
    """Split data into train, test_control, test_stimulated.

    Returns:
        train_adata: All conditions for train cell types
        test_ctrl: Control cells of test cell type (model input)
        test_stim: Stimulated cells of test cell type (ground truth)
    """
    # Training data: all conditions for train cell types
    train_mask = adata.obs[cell_type_key].isin(train_cell_types)
    train_adata = adata[train_mask].copy()

    # Test data: held-out cell type
    test_mask = adata.obs[cell_type_key] == test_cell_type
    test_ctrl = adata[test_mask & (adata.obs[condition_key] == control_label)].copy()
    test_stim = adata[test_mask & (adata.obs[condition_key] == stim_label)].copy()

    print(f"Train: {train_adata.n_obs} cells ({len(train_cell_types)} cell types)")
    print(f"Test control: {test_ctrl.n_obs} cells")
    print(f"Test stimulated (ground truth): {test_stim.n_obs} cells")

    return train_adata, test_ctrl, test_stim


def compute_mean_expression(adata: AnnData) -> np.ndarray:
    """Compute mean gene expression across all cells."""
    if hasattr(adata.X, "toarray"):
        return np.array(adata.X.toarray().mean(axis=0)).flatten()
    return np.array(adata.X.mean(axis=0)).flatten()
