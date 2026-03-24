"""Meta-benchmark orchestrator for all 27 models."""
import scanpy as sc
from pathlib import Path
from typing import List, Dict
import pandas as pd

class MetaBenchmark:
    def __init__(self, adata_path: str, train_cell_types: List[str], test_cell_type: str):
        self.adata = sc.read_h5ad(adata_path)
        self.train_cell_types = train_cell_types
        self.test_cell_type = test_cell_type

    def run_all_models(self, output_dir: str = "results/", n_jobs: int = 1) -> pd.DataFrame:
        """Run all 27 models and return metrics."""
        import sys
        from pathlib import Path as P
        skill_path = P(__file__).parent.parent
        if str(skill_path) not in sys.path:
            sys.path.insert(0, str(skill_path))
        from sub_skills import GenerativeModels, DeepNetModels, ClassicalMLModels, LinearModels

        results = []

        # Module A: Generative (7 models)
        gen_models = GenerativeModels(self.adata)
        results.extend(gen_models.run_all(self.train_cell_types, self.test_cell_type))

        # Module B: Deep Networks (8 models)
        deep_models = DeepNetModels(self.adata)
        results.extend(deep_models.run_all(self.train_cell_types, self.test_cell_type))

        # Module C: Classical ML (6 models)
        ml_models = ClassicalMLModels(self.adata)
        results.extend(ml_models.run_all(self.train_cell_types, self.test_cell_type))

        # Module D: Linear (6 models)
        linear_models = LinearModels(self.adata)
        results.extend(linear_models.run_all(self.train_cell_types, self.test_cell_type))

        df = pd.DataFrame(results)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        df.to_csv(f"{output_dir}/all_models_comparison.csv", index=False)
        return df
