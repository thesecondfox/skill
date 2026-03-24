#!/usr/bin/env python3
"""Example script to run meta-benchmark."""
import sys
sys.path.insert(0, '..')

from scgen_meta_benchmark_skill.core.orchestrator import MetaBenchmark

# Configuration
DATA_PATH = "../data/pbmc_stimulated.h5ad"
TRAIN_CELL_TYPES = ["CD4T", "B_cell", "NK"]
TEST_CELL_TYPE = "CD8T"
OUTPUT_DIR = "../results/"

# Run benchmark
benchmark = MetaBenchmark(
    adata_path=DATA_PATH,
    train_cell_types=TRAIN_CELL_TYPES,
    test_cell_type=TEST_CELL_TYPE
)

results = benchmark.run_all_models(output_dir=OUTPUT_DIR, n_jobs=1)
print(results)
