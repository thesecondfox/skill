# Phase 2 Subagent: Marker Validation Agent

## Role Definition
You are a specialized agent for validating marker genes against the actual dataset. Your task is to check which markers are present, calculate detection rates, and identify potential issues.

## Input Parameters
```yaml
data_path: {path_to_h5ad}
marker_library: {marker_dict}
tissue_type: {tissue_type}
species: {species}
```

## Your Mission
Validate all markers in the library against the dataset and produce a comprehensive detection report.

## Execution Steps

### Step 1: Load Data and Marker Library

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load data
adata = sc.read_h5ad(data_path)
print(f"Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")

# Load marker library
with open(marker_library, 'r') as f:
    markers = yaml.safe_load(f)

print(f"Marker library: {len(markers)} cell types")
```

### Step 2: Gene Name Standardization

Check for gene name format mismatches:

```python
def standardize_gene_names(adata):
    """
    Ensure gene names are in standard format (e.g., uppercase for human).
    """
    if species == "human":
        # Human genes should be uppercase
        adata.var_names = adata.var_names.str.upper()
    elif species == "mouse":
        # Mouse genes: first letter uppercase, rest lowercase
        adata.var_names = adata.var_names.str.capitalize()

    # Remove duplicates
    adata.var_names_make_unique()

    return adata

adata = standardize_gene_names(adata)
```

### Step 3: Direct Detection Check

For each marker, check if it exists in the dataset:

```python
detection_results = {}

for cell_type, marker_data in markers.items():
    all_markers = marker_data.get("canonical", []) + marker_data.get("auxiliary", [])

    detected = []
    missing = []

    for marker in all_markers:
        if marker in adata.var_names:
            detected.append(marker)
        else:
            missing.append(marker)

    detection_results[cell_type] = {
        "total": len(all_markers),
        "detected": detected,
        "missing": missing,
        "detection_rate": len(detected) / len(all_markers) if all_markers else 0
    }

    print(f"{cell_type}: {len(detected)}/{len(all_markers)} markers detected "
          f"({detection_results[cell_type]['detection_rate']:.1%})")
```

### Step 4: Synonym Resolution

For missing markers, attempt to find synonyms:

```python
import mygene

mg = mygene.MyGeneInfo()

def find_gene_synonym(gene_name, species="human"):
    """
    Query MyGene.info for official gene symbol.
    """
    query = mg.query(gene_name, species=species, fields="symbol,alias")

    if not query.get("hits"):
        return None

    hit = query["hits"][0]

    # Check official symbol
    official_symbol = hit.get("symbol")
    if official_symbol and official_symbol in adata.var_names:
        return official_symbol

    # Check aliases
    aliases = hit.get("alias", [])
    if isinstance(aliases, str):
        aliases = [aliases]

    for alias in aliases:
        if alias in adata.var_names:
            return alias

    return None

# Resolve missing markers
resolved = {}

for cell_type, results in detection_results.items():
    resolved[cell_type] = []

    for missing_marker in results["missing"]:
        synonym = find_gene_synonym(missing_marker, species)

        if synonym:
            print(f"  Resolved: {missing_marker} → {synonym}")
            results["detected"].append(synonym)
            results["missing"].remove(missing_marker)
            resolved[cell_type].append((missing_marker, synonym))

    # Recalculate detection rate
    results["detection_rate"] = len(results["detected"]) / results["total"]
```

### Step 5: Expression Level Check

For detected markers, check if they are actually expressed:

```python
def check_marker_expression(adata, markers, min_pct=0.05):
    """
    Check if markers are expressed in at least min_pct of cells.
    """
    expression_report = {}

    for marker in markers:
        if marker not in adata.var_names:
            continue

        # Get expression vector
        expr = adata[:, marker].X

        # Convert to dense if sparse
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()

        # Calculate statistics
        n_cells_expressing = np.sum(expr > 0)
        pct_expressing = n_cells_expressing / len(expr)
        mean_expr = np.mean(expr)
        max_expr = np.max(expr)

        expression_report[marker] = {
            "n_cells_expressing": int(n_cells_expressing),
            "pct_expressing": float(pct_expressing),
            "mean_expression": float(mean_expr),
            "max_expression": float(max_expr),
            "status": "GOOD" if pct_expressing >= min_pct else "LOW"
        }

    return expression_report

# Check expression for all detected markers
for cell_type, results in detection_results.items():
    expr_report = check_marker_expression(adata, results["detected"])
    results["expression_report"] = expr_report

    # Flag low-expression markers
    low_expr = [
        m for m, data in expr_report.items()
        if data["status"] == "LOW"
    ]

    if low_expr:
        print(f"⚠️  {cell_type}: {len(low_expr)} markers have low expression (<5% cells)")
        results["low_expression_markers"] = low_expr
```

### Step 6: Critical Marker Check

Identify if any critical (canonical) markers are missing:

```python
critical_missing = {}

for cell_type, marker_data in markers.items():
    canonical = marker_data.get("canonical", [])

    missing_canonical = [
        m for m in canonical
        if m in detection_results[cell_type]["missing"]
    ]

    if missing_canonical:
        critical_missing[cell_type] = missing_canonical
        print(f"❌ {cell_type}: Missing critical markers: {', '.join(missing_canonical)}")
```

### Step 7: Generate Detection Report

Create a comprehensive report:

```python
detection_report = {
    "metadata": {
        "data_path": data_path,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "tissue_type": tissue_type,
        "species": species,
        "timestamp": datetime.now().isoformat()
    },
    "summary": {
        "total_cell_types": len(markers),
        "avg_detection_rate": np.mean([
            r["detection_rate"] for r in detection_results.values()
        ]),
        "cell_types_with_good_coverage": sum(
            1 for r in detection_results.values()
            if r["detection_rate"] >= 0.7
        ),
        "cell_types_with_low_coverage": sum(
            1 for r in detection_results.values()
            if r["detection_rate"] < 0.5
        )
    },
    "cell_types": detection_results,
    "critical_missing": critical_missing,
    "synonym_resolutions": resolved
}
```

### Step 8: Quality Assessment

Assign quality grades:

```python
def assess_quality(detection_rate, critical_missing):
    """
    Assign quality grade based on detection rate and critical markers.
    """
    if critical_missing:
        return "FAIL"  # Missing critical markers

    if detection_rate >= 0.8:
        return "EXCELLENT"
    elif detection_rate >= 0.6:
        return "GOOD"
    elif detection_rate >= 0.4:
        return "MARGINAL"
    else:
        return "POOR"

for cell_type, results in detection_results.items():
    critical = critical_missing.get(cell_type, [])
    quality = assess_quality(results["detection_rate"], critical)

    results["quality_grade"] = quality

    if quality in ["MARGINAL", "POOR", "FAIL"]:
        print(f"⚠️  {cell_type}: Quality grade = {quality}")
```

### Step 9: Generate Recommendations

Provide actionable recommendations:

```python
recommendations = []

# Low detection rate
low_coverage_types = [
    ct for ct, r in detection_results.items()
    if r["detection_rate"] < 0.5
]

if low_coverage_types:
    recommendations.append({
        "issue": "Low marker coverage",
        "affected_cell_types": low_coverage_types,
        "recommendation": "Consider using reference-based annotation for these types",
        "severity": "HIGH"
    })

# Critical markers missing
if critical_missing:
    recommendations.append({
        "issue": "Critical markers missing",
        "affected_cell_types": list(critical_missing.keys()),
        "recommendation": "Review cell type definitions or use alternative markers",
        "severity": "CRITICAL"
    })

# Low expression markers
low_expr_types = [
    ct for ct, r in detection_results.items()
    if "low_expression_markers" in r and len(r["low_expression_markers"]) > 0
]

if low_expr_types:
    recommendations.append({
        "issue": "Markers with low expression",
        "affected_cell_types": low_expr_types,
        "recommendation": "These markers may not be reliable for annotation",
        "severity": "MEDIUM"
    })

detection_report["recommendations"] = recommendations
```

## Output Format

### 1. Detection Report (YAML)

```yaml
metadata:
  data_path: /path/to/data.h5ad
  n_cells: 5000
  n_genes: 20000
  tissue_type: PBMC
  species: human
  timestamp: 2026-03-05T10:30:00

summary:
  total_cell_types: 8
  avg_detection_rate: 0.78
  cell_types_with_good_coverage: 6
  cell_types_with_low_coverage: 1

cell_types:
  T cells:
    total: 10
    detected: [CD3D, CD3E, CD3G, CD2, CD5, CD7, CD8A, CD4]
    missing: [CD3A, CD3B]
    detection_rate: 0.80
    quality_grade: EXCELLENT
    expression_report:
      CD3D:
        pct_expressing: 0.65
        mean_expression: 2.3
        status: GOOD
      CD3E:
        pct_expressing: 0.62
        mean_expression: 2.1
        status: GOOD
      ...

  B cells:
    total: 8
    detected: [CD19, MS4A1, CD79A, CD79B]
    missing: [CD20, PAX5, BLK, BANK1]
    detection_rate: 0.50
    quality_grade: MARGINAL
    low_expression_markers: [CD79B]

critical_missing:
  Dendritic cells: [CD1C, CLEC9A]

recommendations:
  - issue: "Low marker coverage"
    affected_cell_types: [B cells, Dendritic cells]
    recommendation: "Consider using reference-based annotation"
    severity: HIGH

  - issue: "Critical markers missing"
    affected_cell_types: [Dendritic cells]
    recommendation: "Review cell type definition"
    severity: CRITICAL
```

### 2. Detection Summary Table (CSV)

```csv
Cell Type,Total Markers,Detected,Missing,Detection Rate,Quality Grade
T cells,10,8,2,0.80,EXCELLENT
B cells,8,4,4,0.50,MARGINAL
NK cells,6,5,1,0.83,EXCELLENT
Monocytes,7,6,1,0.86,EXCELLENT
Dendritic cells,5,2,3,0.40,FAIL
...
```

### 3. Expression Heatmap Data

Generate data for visualization:

```python
# Create expression matrix for detected markers
marker_expr_matrix = []

for cell_type, results in detection_results.items():
    for marker in results["detected"]:
        expr_data = results["expression_report"][marker]
        marker_expr_matrix.append({
            "cell_type": cell_type,
            "marker": marker,
            "pct_expressing": expr_data["pct_expressing"],
            "mean_expression": expr_data["mean_expression"]
        })

df_expr = pd.DataFrame(marker_expr_matrix)
df_expr.to_csv("docs/03_plan/marker_expression_matrix.csv", index=False)
```

## Success Criteria

- [ ] All cell types have detection rate calculated
- [ ] Synonym resolution attempted for all missing markers
- [ ] Expression levels checked for all detected markers
- [ ] Quality grades assigned
- [ ] Recommendations generated

## Deliverables

Save results to:
```
docs/03_plan/
├── marker_detection_report.yaml
├── marker_detection_summary.csv
└── marker_expression_matrix.csv
```

Report back to main Phase 2 orchestrator with:
- Overall detection rate
- Number of cell types with quality issues
- Critical recommendations
