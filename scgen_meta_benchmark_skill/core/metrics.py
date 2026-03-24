"""Metrics for evaluating perturbation predictions."""
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score, mean_squared_error
from anndata import AnnData
from typing import Dict, Optional


def compute_metrics(
    predicted: AnnData,
    ground_truth: AnnData,
    gene_names: Optional[list] = None,
    top_n_deg: int = 100,
) -> Dict[str, float]:
    """Compute comprehensive evaluation metrics.

    Args:
        predicted: Predicted stimulated cells
        ground_truth: Actual stimulated cells
        gene_names: Optional gene list for DEG-specific metrics
        top_n_deg: Number of top differentially expressed genes to evaluate

    Returns:
        Dictionary with R2, MSE, Pearson, Spearman, and DEG-specific metrics
    """
    # Extract mean expression vectors
    pred_mean = _get_mean(predicted)
    true_mean = _get_mean(ground_truth)

    # Global metrics
    r2 = r2_score(true_mean, pred_mean)
    mse = mean_squared_error(true_mean, pred_mean)
    pearson_r, pearson_p = pearsonr(true_mean, pred_mean)
    spearman_r, spearman_p = spearmanr(true_mean, pred_mean)

    metrics = {
        "r2_score": float(r2),
        "mse": float(mse),
        "rmse": float(np.sqrt(mse)),
        "pearson_r": float(pearson_r),
        "pearson_pval": float(pearson_p),
        "spearman_r": float(spearman_r),
        "spearman_pval": float(spearman_p),
        "mae": float(np.mean(np.abs(true_mean - pred_mean))),
    }

    # DEG-specific metrics (top N genes by absolute difference in ground truth)
    if top_n_deg > 0:
        diff = np.abs(true_mean - np.mean(true_mean))
        top_idx = np.argsort(diff)[-top_n_deg:]

        deg_r2 = r2_score(true_mean[top_idx], pred_mean[top_idx])
        deg_pearson, _ = pearsonr(true_mean[top_idx], pred_mean[top_idx])

        metrics["deg_r2_top100"] = float(deg_r2)
        metrics["deg_pearson_top100"] = float(deg_pearson)

    return metrics


def compute_per_gene_correlation(predicted: AnnData, ground_truth: AnnData) -> np.ndarray:
    """Compute per-gene Pearson correlation across cells."""
    pred_mat = _to_dense(predicted.X)
    true_mat = _to_dense(ground_truth.X)

    n_genes = pred_mat.shape[1]
    correlations = np.zeros(n_genes)

    for g in range(n_genes):
        if np.std(true_mat[:, g]) > 0 and np.std(pred_mat[:, g]) > 0:
            correlations[g], _ = pearsonr(true_mat[:, g], pred_mat[:, g])
        else:
            correlations[g] = 0.0

    return correlations


def evaluate_prediction(predicted: AnnData, ground_truth: AnnData) -> Dict[str, float]:
    """Convenience wrapper for compute_metrics."""
    return compute_metrics(predicted, ground_truth)


def _get_mean(adata: AnnData) -> np.ndarray:
    """Extract mean expression vector."""
    X = _to_dense(adata.X)
    return np.array(X.mean(axis=0)).flatten()


def _to_dense(X) -> np.ndarray:
    """Convert sparse matrix to dense if needed."""
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)
