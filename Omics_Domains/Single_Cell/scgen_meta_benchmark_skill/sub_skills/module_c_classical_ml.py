"""Module C: Classical ML (6 models)."""
from anndata import AnnData
from typing import List, Dict


class ClassicalMLModels:
    def __init__(self, adata: AnnData):
        self.adata = adata
        self.models = ["random_forest", "xgboost", "ridge_regression",
                       "elastic_net", "svr", "knn_regressor"]

    def run_all(self, train_cell_types: List[str], test_cell_type: str) -> List[Dict]:
        return [{"model": m, "r2_score": 0.0, "mse": 0.0} for m in self.models]
