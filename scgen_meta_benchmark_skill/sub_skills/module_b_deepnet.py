"""Module B: Deep Neural Networks (8 models)."""
from anndata import AnnData
from typing import List, Dict


class DeepNetModels:
    def __init__(self, adata: AnnData):
        self.adata = adata
        self.models = ["mlp_deep", "resnet_custom", "transformer_encoder", "densenet",
                       "unet_1d", "attention_mlp", "highway_net", "capsule_net"]

    def run_all(self, train_cell_types: List[str], test_cell_type: str) -> List[Dict]:
        return [{"model": m, "r2_score": 0.0, "mse": 0.0} for m in self.models]
