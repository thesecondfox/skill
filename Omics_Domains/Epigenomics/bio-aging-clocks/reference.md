# Bio-Aging-Clocks Reference

## Documentation Links

- **GitHub**: https://github.com/bio-learn/biolearn
- **Documentation**: https://bio-learn.github.io/
- **PyPI**: https://pypi.org/project/biolearn/

## Installation

```bash
pip install biolearn
# PyTorch required for AltumAge
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## Core Functions (Sequential Pipeline)

### 1. DataLibrary — Data Ingestion

```python
from biolearn.data_library import DataLibrary, GeoData

data_lib = DataLibrary()

# Browse available sources (2200+)
print(f'Total sources: {len(data_lib.sources)}')
for s in data_lib.sources:
    print(s.id, s.title, s.format, s.organism)

# Filter by organism/format
human_450k = data_lib.lookup_sources(organism='Human', format='Illumina450k')

# Load a GEO dataset (downloads on first call, cached afterward)
source = data_lib.get('GSE30870')   # Returns DataSource
data = source.load()                 # Returns GeoData
```

### 2. GeoData — Data Container

```python
# Constructor signature
GeoData(metadata, dnam=None, rna=None, protein_alamar=None, protein_olink=None)

# Attributes
data.dnam         # pd.DataFrame: CpG sites × samples (beta values 0–1)
data.rna          # pd.DataFrame: genes × samples (expression values)
data.metadata     # pd.DataFrame: sample annotations (age, sex, disease_status, ...)

# Quality assessment
report = data.quality_report()
report.sample_report   # Per-sample quality metrics
```

### 3. Custom Data Loading

```python
# --- Methylation data ---
dnam = pd.read_csv('methylation_matrix.csv', index_col=0)      # CpGs × samples
metadata = pd.read_csv('sample_metadata.csv', index_col=0)     # samples × attributes
custom_data = GeoData(metadata=metadata, dnam=dnam)

# --- From CSV folder (must contain dnam.csv + metadata.csv) ---
custom_data = GeoData.load_csv('/path/to/folder/', name='my_dataset')

# --- RNA-seq data (MUST use rna= parameter) ---
rna_matrix = pd.read_csv('expression_matrix.csv', index_col=0)  # genes × samples
custom_data = GeoData(metadata=metadata, rna=rna_matrix)

# --- Proteomics data ---
custom_data = GeoData(metadata=metadata, protein_olink=protein_df)
```

### 4. ModelGallery — Clock Execution

```python
from biolearn.model_gallery import ModelGallery

gallery = ModelGallery()

# List all models
for name in sorted(gallery.model_definitions.keys()):
    mdef = gallery.model_definitions[name]
    print(name, mdef['model']['type'], mdef.get('species'), mdef.get('tissue'))

# Search models by criteria
blood_models = gallery.search(species='Human', tissue='Blood')

# Run a DNAm clock (standard path)
model = gallery.get('Horvathv1')
results = model.predict(data)
# results: DataFrame with index=sample IDs, column='Predicted'

# Run a transcriptomic clock (MUST bypass ImputationDecorator)
model = gallery.get('Pasta')
results = model.clock.predict(data)  # NOT model.predict(data)
```

### 5. Visualization Functions

```python
from biolearn.visualize import (
    plot_age_prediction,
    plot_clock_correlation_matrix,
    plot_sample_deviations,
    compute_methylation_stats,
    plot_methylation_vs_age,
)

# Predicted vs actual age scatter
models = [gallery.get('Horvathv1'), gallery.get('Hannum'), gallery.get('PhenoAge')]
plot_age_prediction(models, data)

# Clock correlation heatmap
plot_clock_correlation_matrix(models, data)

# Sample quality ridge plot
datasets = {'GSE30870': data}
plot_sample_deviations(datasets)

# CpG-level methylation vs age
combined, stats = compute_methylation_stats(datasets, cpg_site='cg16867657')
plot_methylation_vs_age(combined, stats, cpg_site='cg16867657')
```

## Key Parameters

| Parameter / Field | Context | Description |
|-------------------|---------|-------------|
| `age` | metadata column | Chronological age (numeric). Required for age acceleration. |
| `sex` | metadata column | `"1"` = female, `"2"` = male. Required for GrimAge. |
| `dnam=` | GeoData constructor | Methylation beta-value matrix (CpG × samples). For DNAm clocks. |
| `rna=` | GeoData constructor | Gene expression matrix (genes × samples). For transcriptomic clocks. |
| `protein_olink=` | GeoData constructor | Olink protein matrix. For OrganAge models. |
| `.clock.predict()` | Model execution | Bypass ImputationDecorator for transcriptomic models. |

## Expected Outputs

| Output | Type | Description |
|--------|------|-------------|
| Prediction results | `pd.DataFrame` | Index = sample IDs, column = `'Predicted'` |
| Age acceleration | `pd.Series` | Residual from linear regression of predicted on chronological age |
| Cell fractions | `pd.DataFrame` | Columns = cell types (CD8T, CD4T, NK, Bcell, Mono, Gran, ...) |
| Quality report | `QualityReport` | `.sample_report` attribute with per-sample deviation scores |
| Scatter plot | matplotlib figure | Predicted vs actual age |
| Correlation heatmap | matplotlib figure | Inter-clock Pearson correlations |
