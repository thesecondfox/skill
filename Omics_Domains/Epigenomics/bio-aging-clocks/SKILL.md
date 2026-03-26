---
name: bio-aging-clocks
description: Epigenetic aging clock analysis using biolearn. Use for DNA methylation age prediction (Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE), blood cell deconvolution, disease risk scores, data quality assessment, and biological age visualization. Supports 67+ models across Illumina 27K/450K/EPIC platforms with 2200+ GEO datasets.
license: MIT
metadata:
    skill-author: yzhou
---

# Bio-Aging-Clocks: Epigenetic Age Analysis with biolearn

## Overview

biolearn is a Python library for biological aging analysis, providing 67+ epigenetic clock models and access to 2200+ GEO datasets. It supports DNA methylation-based age prediction, blood cell deconvolution, disease/trait risk scoring, and telomere length estimation. This skill covers the complete workflow from data loading through model prediction and result visualization.

## When to Use This Skill

This skill should be used when:
- Predicting biological/epigenetic age from DNA methylation data
- Running aging clocks (Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE, etc.)
- Loading and analyzing GEO methylation datasets
- Performing blood cell type deconvolution from methylation arrays
- Computing disease/trait risk scores (BMI, smoking, alcohol, CVD, etc.)
- Assessing DNA methylation data quality
- Comparing multiple aging clocks on the same dataset
- Visualizing age acceleration and clock correlations

## Installation

```bash
pip install biolearn
# torch is required for some models (e.g., AltumAge)
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## Quick Start

### Basic Import and Setup

```python
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import DataLibrary
import numpy as np
import pandas as pd
```

### Loading Data from GEO

```python
# DataLibrary provides access to 2200+ curated GEO datasets
data_lib = DataLibrary()

# Get a DataSource object, then call .load() to download and parse
source = data_lib.get('GSE30870')  # Returns DataSource
data = source.load()               # Returns GeoData (downloads on first call)

# GeoData structure
data.dnam       # DataFrame: CpG sites (rows) x samples (columns), shape e.g. (485577, 40)
data.metadata   # DataFrame: sample metadata with columns like 'age', 'sex', 'disease_status'
```

### Running an Aging Clock

```python
gallery = ModelGallery()
model = gallery.get('Horvathv1')
results = model.predict(data)  # Returns DataFrame with 'Predicted' column

# Compare predicted vs chronological age
merged = data.metadata[['age']].join(results, how='inner')
merged.columns = ['Chronological', 'Predicted']
corr = merged['Chronological'].corr(merged['Predicted'])
mae = np.abs(merged['Chronological'] - merged['Predicted']).mean()
print(f'Pearson r = {corr:.4f}, MAE = {mae:.2f} years')
```

## Available Models (67+)

### Epigenetic Age Clocks (DNA Methylation)

| Model | Tissue | Description |
|-------|--------|-------------|
| `Horvathv1` | Multi-tissue | Pan-tissue clock (353 CpGs) |
| `Horvathv2` | Skin + blood | Skin & blood clock |
| `Hannum` | Blood | Blood-based clock (71 CpGs) |
| `PhenoAge` | Blood | Mortality-associated phenotypic age |
| `GrimAgeV1` / `GrimAgeV2` | Blood | Mortality predictor using DNAm surrogates |
| `DunedinPACE` | Blood | Pace of aging (rate, not level) |
| `DunedinPoAm38` | Blood | Earlier pace-of-aging measure |
| `AltumAge` | Multi-tissue | Deep learning clock (PyTorch) |
| `DNAmTL` | Blood, Adipose | Telomere length estimator |
| `EpiTOC1` / `EpiTOC2` | Blood | Mitotic clock (cell division rate) |
| `PCHorvath1` | Multi-tissue | Principal component version of Horvath |
| `HRSInCHPhenoAge` | Blood | PhenoAge variant |
| `DNAmClockCortical` | Human Cortex | Brain-specific clock |
| `Bocklandt` / `Garagnani` / `Weidner` / `Lin` / `VidalBralo` | Blood | Early-generation clocks |
| `Bohlin` / `Knight` | Cord Blood | Gestational age clocks |
| `LeeControl` / `LeeRobust` / `LeeRefinedRobust` / `Mayne` | Placenta | Placental clocks |
| `PEDBE` | Buccal | Pediatric buccal clock |
| `Zhang_10` | Blood | Zhang clock |
| `MiAge` | Blood | Mitotic age |
| `GPAge10/30/71/A/B/C` | Blood | Gaussian process age models |
| `YingAdaptAge` / `YingCausAge` / `YingDamAge` | Blood | Ying decomposed aging clocks |
| `StocH` / `StocP` / `StocZ` | Multi-tissue/Blood | Stochastic epigenetic clocks |

### Disease & Trait Risk Scores

| Model | What it predicts |
|-------|-----------------|
| `BMI_McCartney` / `BMI_Reed` | Body mass index |
| `SmokingMcCartney` | Smoking exposure |
| `AlcoholMcCartney` | Alcohol consumption |
| `CVD_Westerman` | Cardiovascular disease risk |
| `AD_Bahado-Singh` | Alzheimer's disease risk |
| `DepressionBarbu` | Depression risk |
| `DownSyndrome` | Down syndrome |
| `ProstateCancerKirby` | Prostate cancer risk |
| `HepatoXu` | Liver (circulating DNA) |
| `BodyFatMcCartney` | Body fat percentage |
| `EducationMcCartney` | Educational attainment |
| `HDLCholesterolMcCartney` / `LDLCholesterolMcCartney` / `TotalCholesterolMcCartney` | Cholesterol levels |

### Blood Cell Deconvolution

| Model | Platform |
|-------|----------|
| `DeconvoluteBlood450K` | Illumina 450K |
| `DeconvoluteBloodEPIC` | Illumina EPIC |
| `TwelveCellDeconvoluteBloodEPIC` | EPIC (12 cell types) |

### Transcriptomic & Proteomic Clocks

| Model | Type |
|-------|------|
| `Pasta` / `REG` / `TranscriptomicPredictionModel` | RNA-based age |
| `OrganAgeChronological` / `OrganAgeMortality` | Proteomic organ age |
| `OrganAge1500Chronological` / `OrganAge1500Mortality` | Proteomic organ age (1500 features) |

### Other

| Model | Description |
|-------|-------------|
| `SexEstimation` | Predict sex from methylation (output columns: X, Y, predicted_sex) |
| `HurdleInflammAge` | Inflammation age (API-based) |

## Core API Reference

### DataLibrary

```python
from biolearn.data_library import DataLibrary, GeoData

data_lib = DataLibrary()

# Browse available sources
print(f'Total sources: {len(data_lib.sources)}')
for s in data_lib.sources:
    print(s.id, s.title, s.format, s.organism)

# Filter by organism/format
human_450k = data_lib.lookup_sources(organism='Human', format='Illumina450k')

# Load a dataset
source = data_lib.get('GSE40279')  # Returns DataSource
data = source.load()                # Returns GeoData (may take time for large datasets)
```

### GeoData Structure

```python
data.dnam         # pd.DataFrame: CpG sites x samples (methylation beta values 0-1)
data.metadata     # pd.DataFrame: sample annotations (age, sex, disease_status, etc.)

# Quality assessment
report = data.quality_report()
report.sample_report   # Per-sample quality metrics including deviation scores
```

### ModelGallery

```python
from biolearn.model_gallery import ModelGallery

gallery = ModelGallery()

# List all models
for name in sorted(gallery.model_definitions.keys()):
    mdef = gallery.model_definitions[name]
    print(name, mdef['model']['type'], mdef.get('species'), mdef.get('tissue'))

# Search models by criteria
blood_models = gallery.search(species='Human', tissue='Blood')

# Get and run a model
model = gallery.get('Horvathv1')
results = model.predict(data)  # data is a GeoData object
# results: DataFrame with index=sample IDs, column='Predicted'
```

## Clock Selection Guide

Choosing the right clock depends on your research question, tissue type, and data platform. Below is a decision framework.

### By Research Question

| Goal | Recommended Clocks | Rationale |
|------|-------------------|-----------|
| Chronological age prediction | Horvathv1, Hannum, PCHorvath1 | First-generation clocks trained directly on chronological age |
| Biological age / mortality risk | PhenoAge, GrimAgeV1/V2 | Trained on mortality-associated biomarkers, better at capturing health-relevant aging |
| Pace of aging (rate) | DunedinPACE, DunedinPoAm38 | Measures how fast someone is aging (rate of change), not biological age level |
| Cellular replication / mitotic age | EpiTOC1, EpiTOC2, MiAge | Estimates cumulative cell divisions, relevant for cancer risk |
| Telomere length | DNAmTL | Epigenetic surrogate for leukocyte telomere length |
| Gestational age | Bohlin, Knight | Designed for cord blood samples |
| Placental age | LeeControl, LeeRobust, LeeRefinedRobust, Mayne | Placenta-specific clocks |
| Pediatric age | PEDBE | Buccal epithelial cells in children |
| Brain aging | DNAmClockCortical | Human cortex tissue only |
| Organ-specific aging | OrganAgeChronological, OrganAgeMortality | Proteomic (Olink) organ age clocks |
| Comprehensive profiling | Run multiple clocks | Compare 1st-gen (Horvath/Hannum) + 2nd-gen (PhenoAge/GrimAge) + pace (DunedinPACE) |

### By Tissue Type

| Tissue | Compatible Clocks |
|--------|------------------|
| Blood (most common) | Horvathv1, Hannum, PhenoAge, GrimAgeV1/V2, DunedinPACE, DNAmTL, EpiTOC1/2, GPAge series, all McCartney scores |
| Multi-tissue | Horvathv1, Horvathv2, AltumAge, StocH, PCHorvath1 |
| Skin + Blood | Horvathv2 |
| Cord Blood | Bohlin, Knight |
| Placenta | LeeControl, LeeRobust, LeeRefinedRobust, Mayne |
| Buccal | PEDBE |
| Human Cortex | DNAmClockCortical |
| Circulating DNA | HepatoXu |
| Prostate | ProstateCancerKirby |

### By Platform

| Platform | Notes |
|----------|-------|
| Illumina 450K | Most clocks work; use DeconvoluteBlood450K for deconvolution |
| Illumina EPIC | Most clocks work; use DeconvoluteBloodEPIC or TwelveCellDeconvoluteBloodEPIC |
| Illumina 27K | Only clocks using ≤27K CpGs (e.g., Horvathv1 partially, Bocklandt, Garagnani, Weidner) |
| RNA-seq | Pasta, REG, TranscriptomicPredictionModel only |
| Olink proteomics | OrganAge models only |

### By Generation / Design Philosophy

| Generation | Clocks | Characteristics |
|-----------|--------|----------------|
| 1st gen (age-trained) | Horvathv1, Hannum, Bocklandt, Garagnani, Weidner, Lin | Trained on chronological age; good correlation but less health-predictive |
| 2nd gen (outcome-trained) | PhenoAge, GrimAgeV1/V2 | Trained on mortality/health biomarkers; better at predicting health outcomes |
| 3rd gen (pace) | DunedinPACE, DunedinPoAm38 | Measures rate of aging from longitudinal data; complementary to level-based clocks |
| Deep learning | AltumAge | Neural network; requires PyTorch; multi-tissue |
| PC-based | PCHorvath1 | Principal component transformation; more robust to noise |
| Stochastic | StocH, StocP, StocZ | Capture stochastic epigenetic variation |
| Decomposed | YingAdaptAge, YingCausAge, YingDamAge | Separate adaptive, causal, and damage components of aging |

### Recommended Combinations for Common Study Designs

```python
# Aging cohort study — comprehensive panel
clocks_comprehensive = [
    'Horvathv1',      # 1st gen reference
    'Hannum',         # 1st gen blood
    'PhenoAge',       # 2nd gen mortality
    'GrimAgeV1',      # 2nd gen composite
    'DunedinPACE',    # 3rd gen pace
    'DNAmTL',         # telomere length
    'EpiTOC2',        # mitotic age
]

# Quick screening — minimal but informative
clocks_minimal = ['Horvathv1', 'PhenoAge', 'DunedinPACE']

# Cancer-focused
clocks_cancer = ['EpiTOC1', 'EpiTOC2', 'MiAge', 'PhenoAge']

# Cardiovascular / metabolic
clocks_cvd = ['GrimAgeV1', 'PhenoAge', 'CVD_Westerman', 'BMI_McCartney']

# Lifestyle exposure
clocks_lifestyle = ['SmokingMcCartney', 'AlcoholMcCartney', 'BMI_McCartney', 'BodyFatMcCartney', 'EducationMcCartney']

# Prenatal / pediatric
clocks_prenatal = ['Bohlin', 'Knight', 'PEDBE']  # cord blood / buccal
clocks_placental = ['LeeControl', 'LeeRobust', 'LeeRefinedRobust', 'Mayne']
```

### Special Considerations

- **GrimAge** requires `sex` in metadata (standardized as "1"=female, "2"=male)
- **DunedinPACE** measures pace (rate), not level — it does not predict age, but how fast aging is occurring; mean ~1.0 in healthy populations
- **AltumAge** requires PyTorch installed
- **HurdleInflammAge** calls an external API — requires internet access
- **OrganAge** models need Olink proteomic data, not methylation
- **Transcriptomic models** (Pasta, REG) need RNA-seq expression, not methylation
- When comparing clocks, always report which generation/type each belongs to, as they measure fundamentally different aspects of aging

## Common Workflows

### 1. Multi-Clock Comparison

```python
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import DataLibrary
import pandas as pd
import numpy as np

data = DataLibrary().get('GSE30870').load()
gallery = ModelGallery()

clock_names = ['Horvathv1', 'Hannum', 'PhenoAge', 'DNAmTL']
all_predictions = data.metadata[['age']].copy()

for name in clock_names:
    model = gallery.get(name)
    pred = model.predict(data)
    all_predictions[name] = pred['Predicted']

# Compute age acceleration (residual from regression)
from sklearn.linear_model import LinearRegression
for name in clock_names:
    valid = all_predictions[['age', name]].dropna()
    lr = LinearRegression().fit(valid[['age']], valid[name])
    all_predictions[f'{name}_accel'] = all_predictions[name] - lr.predict(all_predictions[['age']])

print(all_predictions.describe())
```

### 2. Blood Cell Deconvolution

```python
model = gallery.get('DeconvoluteBloodEPIC')  # or DeconvoluteBlood450K
cell_fractions = model.predict(data)
print(cell_fractions.head())
# Columns: CD8T, CD4T, NK, Bcell, Mono, Gran (or 12 types for TwelveCellDeconvoluteBloodEPIC)
```

### 3. Data Quality Assessment

```python
from biolearn.visualize import plot_sample_deviations

datasets = {'GSE30870': data}
plot_sample_deviations(datasets)  # Ridge density plot of sample deviations
```

### 4. Age Prediction Visualization

```python
from biolearn.visualize import plot_age_prediction, plot_clock_correlation_matrix

models = [gallery.get('Horvathv1'), gallery.get('Hannum'), gallery.get('PhenoAge')]
plot_age_prediction(models, data)           # Predicted vs actual age scatter
plot_clock_correlation_matrix(models, data) # Clock correlation heatmap
```

### 5. CpG Site Analysis

```python
from biolearn.visualize import compute_methylation_stats, plot_methylation_vs_age

datasets = {'GSE30870': data}
combined, stats = compute_methylation_stats(datasets, cpg_site='cg16867657')
plot_methylation_vs_age(combined, stats, cpg_site='cg16867657')
```

### 6. Loading Custom Data

```python
from biolearn.data_library import GeoData

# From CSV files in a folder (must have dnam.csv and metadata.csv)
custom_data = GeoData.load_csv('/path/to/folder/', name='my_dataset')

# Or construct manually — methylation data
import pandas as pd
dnam = pd.read_csv('methylation_matrix.csv', index_col=0)      # CpGs x samples
metadata = pd.read_csv('sample_metadata.csv', index_col=0)     # samples x attributes
custom_data = GeoData(metadata=metadata, dnam=dnam)

# Run any model on custom data
results = gallery.get('PhenoAge').predict(custom_data)
```

**IMPORTANT: For transcriptomic data, use the `rna=` parameter, NOT `dnam=`:**

```python
# Construct GeoData with RNA-seq expression data
rna_matrix = pd.read_csv('expression_matrix.csv', index_col=0)  # genes x samples
custom_data = GeoData(metadata=metadata, rna=rna_matrix)
# GeoData constructor signature: GeoData(metadata, dnam=None, rna=None, protein_alamar=None, protein_olink=None)
```

### 7. Running Transcriptomic Clocks (Pasta, REG, TranscriptomicPredictionModel)

Transcriptomic clocks have several known issues that require workarounds:

```python
import os
import pandas as pd
import numpy as np
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData

# --- Load RNA-seq data (e.g., TCGA) ---
counts = pd.read_csv('TCGA-COAD.star_counts.tsv.gz', sep='\t', compression='gzip')

# Strip Ensembl version numbers: ENSG00000000003.15 -> ENSG00000000003
# (Pasta/REG coefficients use versionless Ensembl IDs)
counts['Ensembl_ID'] = counts['Ensembl_ID'].str.split('.').str[0]
counts = counts.set_index('Ensembl_ID').groupby(level=0).mean()

metadata = pd.read_csv('clinical.tsv.gz', sep='\t', compression='gzip')
# ... prepare metadata with 'age' column, indexed by sample ID ...

# CRITICAL: use rna= not dnam=
geo_data = GeoData(metadata=metadata, rna=counts)

# --- Workaround 1: Bypass ImputationDecorator ---
# gallery.get() returns an ImputationDecorator that calls .methylation_sites()
# which does NOT exist on LinearTranscriptomicModel → raises AttributeError.
# Fix: call .clock.predict() to bypass the decorator.
gallery = ModelGallery()
model = gallery.get('Pasta')
results = model.clock.predict(geo_data)  # NOT model.predict(geo_data)

# --- Workaround 2: Pasta/REG relative path bug ---
# Pasta's preprocess function reads "biolearn/data/Pasta.csv" using a RELATIVE path.
# This fails unless cwd is the site-packages directory.
# Fix: temporarily chdir before calling predict.
import biolearn
biolearn_parent = os.path.dirname(os.path.dirname(biolearn.__file__))
original_cwd = os.getcwd()
os.chdir(biolearn_parent)  # e.g., .../lib/python3.10/site-packages

for name in ['Pasta', 'REG', 'TranscriptomicPredictionModel']:
    pred = gallery.get(name).clock.predict(geo_data)
    print(name, pred['Predicted'].describe())

os.chdir(original_cwd)
```

**Transcriptomic clock caveats:**
- **Pasta** is sensitive to expression value scale. If the model expects a different normalization (e.g., TPM, raw counts) than your input data, predictions can be wildly off (negative values, MAE > 100 years). Verify your data's normalization matches the model's expectation.

## Key Parameters and Notes

### Dataset Size Considerations
- Large datasets (e.g., GSE40279 with 656 samples x 485K CpGs) can take significant time and memory to download and parse
- First load downloads from GEO; subsequent loads use local cache
- Consider starting with smaller datasets for testing

### Model Compatibility
- Most clocks require Illumina 450K or EPIC methylation data
- Models handle missing CpG sites via imputation internally
- Transcriptomic models (Pasta, REG) need RNA-seq expression data
- Proteomic models (OrganAge) need Olink protein data
- GrimAge models require sex information in metadata

### Metadata Requirements
- `age`: chronological age (numeric) — needed for age acceleration calculation
- `sex`: biological sex (standardized to "1"=female, "2"=male) — needed for GrimAge, SexEstimation
- Additional columns vary by dataset (disease_status, tissue, cell_type, etc.)

## Common Pitfalls

1. **DataSource vs GeoData**: `DataLibrary().get()` returns `DataSource`; call `.load()` to get `GeoData`
2. **Large dataset downloads**: First load of large GEO datasets can take 10+ minutes
3. **Missing torch**: AltumAge and some models require PyTorch — install with `pip install torch`
4. **Memory usage**: 450K/EPIC arrays with many samples consume significant RAM
5. **Platform mismatch**: Ensure your data platform matches the model's expected platform
6. **Sex encoding**: biolearn standardizes sex as "1" (female) and "2" (male)
7. **GeoData `rna=` vs `dnam=`**: Transcriptomic models read from `geo_data.rna`, NOT `geo_data.dnam`. If you pass expression data via `dnam=`, transcriptomic clocks will silently get `None` and crash. Always use `GeoData(metadata=..., rna=expression_df)` for RNA-seq data.
8. **ImputationDecorator breaks transcriptomic models**: `ModelGallery().get()` wraps ALL models in `ImputationDecorator`, which calls `.methylation_sites()` — a method that does not exist on `LinearTranscriptomicModel`. Workaround: use `model.clock.predict(data)` instead of `model.predict(data)` for Pasta, REG, and TranscriptomicPredictionModel.
9. **Pasta/REG relative path bug**: The `preprocess_pasta` function reads `"biolearn/data/Pasta.csv"` as a relative path. This fails unless your working directory is the `site-packages` directory. Workaround: `os.chdir()` to `os.path.dirname(os.path.dirname(biolearn.__file__))` before calling predict, then restore cwd afterward.
10. **Ensembl ID version stripping**: Transcriptomic clock coefficients use versionless Ensembl IDs (e.g., `ENSG00000196839`). If your data has versioned IDs (e.g., `ENSG00000196839.12`), strip versions with `df.index = df.index.str.split('.').str[0]` before prediction.

## Additional Resources

- **GitHub**: https://github.com/bio-learn/biolearn
- **Documentation**: https://bio-learn.github.io/
- **PyPI**: https://pypi.org/project/biolearn/
- **Supported GEO datasets**: Use `DataLibrary().sources` to browse all 2200+ datasets
