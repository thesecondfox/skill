#!/bin/bash
# Setup and run benchmark on cluster

# Load modules
module load python/3.9
module load cuda/11.8

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install scanpy scgen scvi-tools scikit-learn xgboost optuna

# Run benchmark
python run_script.py
