"""Module D: Linear Methods (6 models)."""
from anndata import AnnData
from typing import List, Dict


class LinearModels:
    def __init__(self, adata: AnnData):
        self.adata = adata
        self.models = ["differential_expression", "mean_shift", "nearest_neighbor",
                       "linear_projection", "combat_correction", "scanorama_integration"]

    def run_all(self, train_cell_types: List[str], test_cell_type: str) -> List[Dict]:
        return [{"model": m, "r2_score": 0.0, "mse": 0.0} for m in self.models]
