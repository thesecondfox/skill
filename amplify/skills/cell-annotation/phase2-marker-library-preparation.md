# Phase 2: Marker Gene Library Preparation

## Overview
This phase builds a comprehensive, tissue-specific marker gene library that will be used for annotation. The library combines curated databases, literature sources, and user-provided markers into a unified, prioritized resource.

## Prerequisites
- Phase 1 completed with approved annotation strategy
- `annotation-strategy.yaml` exists with expected cell types
- Data object loaded with gene expression matrix

## Execution Flow

```
Step 1: Automatic Database Retrieval
  ├─ Query CellMarker database
  ├─ Query PanglaoDB
  └─ Query additional tissue-specific databases

Step 2: User-Provided Marker Integration
  ├─ Load custom marker file (if provided)
  ├─ Validate gene names
  └─ Merge with database markers

Step 3: Marker Detection Validation
  ├─ Check which markers are present in dataset
  ├─ Calculate detection rates
  └─ Flag missing critical markers

Step 4: Marker Prioritization
  ├─ Classify: canonical vs auxiliary
  ├─ Score by specificity
  └─ Rank by reliability

Step 5: Conflict Resolution
  ├─ Identify ambiguous markers (expressed in multiple types)
  ├─ Create marker combinations for disambiguation
  └─ Document marker logic

Step 6: Library Export
  └─ Generate `marker_gene_library.yaml`
```

## Step 1: Automatic Database Retrieval

### 1.1 CellMarker Database Query

**Data Source**: http://xteam.xbio.top/CellMarker/

```python
def query_cellmarker(tissue_type, species="human"):
    """
    Query CellMarker database for tissue-specific markers.

    Parameters:
    -----------
    tissue_type : str
        Tissue name (e.g., "PBMC", "brain", "liver")
    species : str
        "human" or "mouse"

    Returns:
    --------
    dict : {cell_type: [marker_genes]}
    """
    import requests
    import pandas as pd

    # CellMarker API endpoint (example - adjust to actual API)
    url = "http://xteam.xbio.top/CellMarker/api/query"

    params = {
        "tissue": tissue_type,
        "species": species,
        "format": "json"
    }

    response = requests.get(url, params=params)

    if response.status_code != 200:
        print(f"Warning: CellMarker query failed with status {response.status_code}")
        return {}

    data = response.json()

    # Parse response into structured format
    markers = {}
    for entry in data.get("results", []):
        cell_type = entry["cell_type"]
        marker_list = entry["markers"].split(",")

        if cell_type not in markers:
            markers[cell_type] = []

        markers[cell_type].extend([m.strip() for m in marker_list])

    # Deduplicate
    markers = {k: list(set(v)) for k, v in markers.items()}

    return markers
```

### 1.2 PanglaoDB Query

**Data Source**: https://panglaodb.se/

```python
def query_panglaodb(tissue_type, species="human"):
    """
    Query PanglaoDB for cell type markers.

    PanglaoDB provides a downloadable marker table.
    Download from: https://panglaodb.se/markers.html
    """
    import pandas as pd

    # Load pre-downloaded PanglaoDB marker table
    # (In practice, download this file once and cache it)
    panglaodb_file = "~/.cache/amplify/panglaodb_markers.tsv"

    try:
        df = pd.read_csv(panglaodb_file, sep="\t")
    except FileNotFoundError:
        print(f"Warning: PanglaoDB file not found at {panglaodb_file}")
        print("Download from: https://panglaodb.se/markers.html")
        return {}

    # Filter by species and tissue
    species_map = {"human": "Hs", "mouse": "Mm"}
    df_filtered = df[
        (df["species"] == species_map[species]) &
        (df["organ"].str.contains(tissue_type, case=False, na=False))
    ]

    # Group by cell type
    markers = {}
    for cell_type, group in df_filtered.groupby("cell type"):
        markers[cell_type] = group["official gene symbol"].tolist()

    return markers
```

### 1.3 Tissue-Specific Database Dispatch

For certain tissues, use specialized databases:

```python
TISSUE_SPECIFIC_DATABASES = {
    "brain": ["BrainCellData", "Allen Brain Atlas"],
    "immune": ["ImmGen", "DICE"],
    "tumor": ["CancerSEA", "TISCH"],
    "liver": ["LiverAtlas"],
}

def query_specialized_db(tissue_type, species="human"):
    """
    Query tissue-specific databases if available.
    """
    if tissue_type not in TISSUE_SPECIFIC_DATABASES:
        return {}

    markers = {}
    for db_name in TISSUE_SPECIFIC_DATABASES[tissue_type]:
        if db_name == "ImmGen":
            markers.update(query_immgen(species))
        elif db_name == "Allen Brain Atlas":
            markers.update(query_allen_brain(species))
        # Add more as needed

    return markers
```

### 1.4 Consolidate Database Results

```python
def consolidate_database_markers(tissue_type, species="human"):
    """
    Query all databases and merge results.
    """
    all_markers = {}

    # Query each database
    sources = {
        "CellMarker": query_cellmarker(tissue_type, species),
        "PanglaoDB": query_panglaodb(tissue_type, species),
        "Specialized": query_specialized_db(tissue_type, species),
    }

    # Merge markers by cell type
    for source_name, markers in sources.items():
        for cell_type, marker_list in markers.items():
            if cell_type not in all_markers:
                all_markers[cell_type] = {
                    "markers": [],
                    "sources": []
                }

            for marker in marker_list:
                if marker not in all_markers[cell_type]["markers"]:
                    all_markers[cell_type]["markers"].append(marker)
                    all_markers[cell_type]["sources"].append(source_name)
                else:
                    # Marker found in multiple sources — increase confidence
                    idx = all_markers[cell_type]["markers"].index(marker)
                    all_markers[cell_type]["sources"][idx] += f", {source_name}"

    return all_markers
```

## Step 2: User-Provided Marker Integration

### 2.1 Load Custom Markers

```python
def load_custom_markers(custom_marker_file):
    """
    Load user-provided marker genes.

    Expected format (YAML):
    ---
    T cells:
      canonical: [CD3D, CD3E, CD3G]
      auxiliary: [CD2, CD5, CD7]
    B cells:
      canonical: [CD19, MS4A1, CD79A]
      auxiliary: [CD79B, BLK]
    """
    import yaml

    if not custom_marker_file:
        return {}

    with open(custom_marker_file, 'r') as f:
        custom_markers = yaml.safe_load(f)

    # Validate structure
    for cell_type, markers in custom_markers.items():
        if not isinstance(markers, dict):
            raise ValueError(
                f"Invalid format for {cell_type}. "
                "Expected dict with 'canonical' and 'auxiliary' keys."
            )

        if "canonical" not in markers:
            raise ValueError(f"Missing 'canonical' markers for {cell_type}")

    return custom_markers
```

### 2.2 Gene Name Validation

```python
def validate_gene_names(marker_dict, adata):
    """
    Check if marker genes exist in the dataset.
    Convert synonyms if needed.
    """
    import mygene

    mg = mygene.MyGeneInfo()
    dataset_genes = set(adata.var_names)

    validated_markers = {}
    missing_markers = {}

    for cell_type, markers in marker_dict.items():
        validated_markers[cell_type] = []
        missing_markers[cell_type] = []

        for marker in markers:
            if marker in dataset_genes:
                validated_markers[cell_type].append(marker)
            else:
                # Try to find synonym
                query_result = mg.query(marker, species="human", fields="symbol")

                if query_result.get("hits"):
                    official_symbol = query_result["hits"][0].get("symbol")
                    if official_symbol in dataset_genes:
                        validated_markers[cell_type].append(official_symbol)
                        print(f"Converted {marker} → {official_symbol}")
                    else:
                        missing_markers[cell_type].append(marker)
                else:
                    missing_markers[cell_type].append(marker)

    return validated_markers, missing_markers
```

### 2.3 Merge Database and Custom Markers

```python
def merge_marker_sources(database_markers, custom_markers):
    """
    Merge database and user-provided markers.
    User markers take priority (marked as 'canonical').
    """
    merged = {}

    # Start with database markers
    for cell_type, data in database_markers.items():
        merged[cell_type] = {
            "canonical": [],
            "auxiliary": data["markers"].copy(),
            "sources": data["sources"].copy()
        }

    # Override with custom markers
    for cell_type, markers in custom_markers.items():
        if cell_type not in merged:
            merged[cell_type] = {
                "canonical": [],
                "auxiliary": [],
                "sources": []
            }

        # User-specified canonical markers
        merged[cell_type]["canonical"] = markers.get("canonical", [])

        # User-specified auxiliary markers
        if "auxiliary" in markers:
            merged[cell_type]["auxiliary"].extend(markers["auxiliary"])

        # Mark source as user-provided
        merged[cell_type]["sources"].append("user_provided")

    # Deduplicate
    for cell_type in merged:
        merged[cell_type]["canonical"] = list(set(merged[cell_type]["canonical"]))
        merged[cell_type]["auxiliary"] = list(set(merged[cell_type]["auxiliary"]))

        # Remove canonical markers from auxiliary
        merged[cell_type]["auxiliary"] = [
            m for m in merged[cell_type]["auxiliary"]
            if m not in merged[cell_type]["canonical"]
        ]

    return merged
```

## Step 3: Marker Detection Validation

### 3.1 Calculate Detection Rates

```python
def calculate_marker_detection(adata, marker_library):
    """
    For each cell type, calculate what % of markers are detected.
    """
    detection_report = {}

    for cell_type, markers in marker_library.items():
        all_markers = markers["canonical"] + markers["auxiliary"]

        detected = [m for m in all_markers if m in adata.var_names]
        missing = [m for m in all_markers if m not in adata.var_names]

        detection_rate = len(detected) / len(all_markers) if all_markers else 0

        detection_report[cell_type] = {
            "total_markers": len(all_markers),
            "detected": len(detected),
            "missing": len(missing),
            "detection_rate": detection_rate,
            "missing_markers": missing,
            "detected_markers": detected
        }

    return detection_report
```

### 3.2 Flag Critical Missing Markers

```python
def flag_critical_missing(detection_report, threshold=0.5):
    """
    Identify cell types with insufficient marker coverage.
    """
    warnings = []

    for cell_type, report in detection_report.items():
        if report["detection_rate"] < threshold:
            warnings.append({
                "cell_type": cell_type,
                "detection_rate": report["detection_rate"],
                "missing_markers": report["missing_markers"],
                "severity": "HIGH" if report["detection_rate"] < 0.3 else "MEDIUM"
            })

    return warnings
```

## Step 4: Marker Prioritization

### 4.1 Specificity Scoring

```python
def score_marker_specificity(adata, marker_library):
    """
    Calculate how specific each marker is to its cell type.

    Specificity = (expression in target) / (expression in all cells)
    """
    import scanpy as sc
    import numpy as np

    specificity_scores = {}

    for cell_type, markers in marker_library.items():
        specificity_scores[cell_type] = {}

        for marker in markers["canonical"] + markers["auxiliary"]:
            if marker not in adata.var_names:
                continue

            # Get expression vector
            expr = adata[:, marker].X.toarray().flatten()

            # Calculate specificity (simplified version)
            # In practice, this would compare expression across clusters
            mean_expr = np.mean(expr)
            max_expr = np.max(expr)

            specificity = mean_expr / (max_expr + 1e-6)

            specificity_scores[cell_type][marker] = {
                "specificity": specificity,
                "mean_expression": mean_expr,
                "max_expression": max_expr
            }

    return specificity_scores
```

### 4.2 Rank Markers by Reliability

```python
def rank_markers(marker_library, specificity_scores, detection_report):
    """
    Rank markers by:
    1. Detection (present in dataset)
    2. Specificity (high = better)
    3. Source reliability (multiple sources = better)
    """
    ranked_markers = {}

    for cell_type, markers in marker_library.items():
        ranked_markers[cell_type] = {
            "canonical": [],
            "auxiliary": []
        }

        for category in ["canonical", "auxiliary"]:
            marker_scores = []

            for marker in markers[category]:
                # Skip if not detected
                if marker not in detection_report[cell_type]["detected_markers"]:
                    continue

                # Get specificity
                spec = specificity_scores[cell_type].get(marker, {}).get("specificity", 0)

                # Count sources
                source_count = markers["sources"].count(marker)

                # Composite score
                score = spec * (1 + 0.1 * source_count)

                marker_scores.append((marker, score))

            # Sort by score (descending)
            marker_scores.sort(key=lambda x: x[1], reverse=True)

            ranked_markers[cell_type][category] = [m for m, s in marker_scores]

    return ranked_markers
```

## Step 5: Conflict Resolution

### 5.1 Identify Ambiguous Markers

```python
def identify_ambiguous_markers(marker_library):
    """
    Find markers that appear in multiple cell types.
    """
    marker_to_celltypes = {}

    for cell_type, markers in marker_library.items():
        for marker in markers["canonical"] + markers["auxiliary"]:
            if marker not in marker_to_celltypes:
                marker_to_celltypes[marker] = []
            marker_to_celltypes[marker].append(cell_type)

    # Filter to ambiguous markers (appear in 2+ types)
    ambiguous = {
        marker: celltypes
        for marker, celltypes in marker_to_celltypes.items()
        if len(celltypes) > 1
    }

    return ambiguous
```

### 5.2 Create Marker Combinations

```python
def create_marker_combinations(marker_library, ambiguous_markers):
    """
    For ambiguous markers, create combinations for disambiguation.

    Example:
    CD4 is expressed in both T cells and monocytes.
    Combination: CD4+ AND CD3+ → T cells
                 CD4+ AND CD14+ → monocytes
    """
    combinations = {}

    for cell_type, markers in marker_library.items():
        canonical = markers["canonical"]

        if len(canonical) < 2:
            continue  # Need at least 2 markers for combination

        # Create pairwise combinations
        from itertools import combinations as comb

        for pair in comb(canonical, 2):
            combo_key = f"{pair[0]}+ AND {pair[1]}+"
            combinations[combo_key] = cell_type

    return combinations
```

### 5.3 Document Marker Logic

```python
def document_marker_logic(marker_library, ambiguous_markers, combinations):
    """
    Create a human-readable logic document.
    """
    logic_doc = []

    logic_doc.append("# Marker Gene Logic\n")
    logic_doc.append("## Unambiguous Markers\n")

    for cell_type, markers in marker_library.items():
        unambiguous = [
            m for m in markers["canonical"]
            if m not in ambiguous_markers
        ]

        if unambiguous:
            logic_doc.append(f"### {cell_type}")
            logic_doc.append(f"- Specific markers: {', '.join(unambiguous)}")

    logic_doc.append("\n## Ambiguous Markers (require combinations)\n")

    for marker, celltypes in ambiguous_markers.items():
        logic_doc.append(f"### {marker}")
        logic_doc.append(f"- Found in: {', '.join(celltypes)}")
        logic_doc.append("- Disambiguation:")

        for combo, celltype in combinations.items():
            if marker in combo:
                logic_doc.append(f"  - {combo} → {celltype}")

    return "\n".join(logic_doc)
```

## Step 6: Library Export

### 6.1 Generate Final Library File

```python
def export_marker_library(
    marker_library,
    ranked_markers,
    detection_report,
    ambiguous_markers,
    combinations,
    output_path="docs/03_plan/marker_gene_library.yaml"
):
    """
    Export the final marker library in YAML format.
    """
    import yaml

    library = {
        "metadata": {
            "created_at": datetime.now().isoformat(),
            "tissue_type": tissue_type,
            "species": species,
            "total_cell_types": len(marker_library),
        },
        "cell_types": {}
    }

    for cell_type, markers in marker_library.items():
        library["cell_types"][cell_type] = {
            "canonical_markers": ranked_markers[cell_type]["canonical"],
            "auxiliary_markers": ranked_markers[cell_type]["auxiliary"],
            "detection_rate": detection_report[cell_type]["detection_rate"],
            "missing_markers": detection_report[cell_type]["missing_markers"],
            "sources": list(set(markers["sources"])),
            "ambiguous_markers": [
                m for m in markers["canonical"]
                if m in ambiguous_markers
            ],
            "marker_combinations": [
                combo for combo, ct in combinations.items()
                if ct == cell_type
            ]
        }

    # Write to file
    with open(output_path, 'w') as f:
        yaml.dump(library, f, default_flow_style=False, sort_keys=False)

    print(f"Marker library exported to: {output_path}")

    return library
```

### 6.2 Generate Summary Report

```python
def generate_summary_report(marker_library, detection_report, warnings):
    """
    Create a human-readable summary for the user.
    """
    report = []

    report.append("# Phase 2 Complete: Marker Gene Library Prepared\n")

    report.append("## Summary Statistics")
    report.append(f"- Total cell types: {len(marker_library)}")

    total_markers = sum(
        len(m["canonical"]) + len(m["auxiliary"])
        for m in marker_library.values()
    )
    report.append(f"- Total markers: {total_markers}")

    avg_detection = sum(
        r["detection_rate"] for r in detection_report.values()
    ) / len(detection_report)
    report.append(f"- Average detection rate: {avg_detection:.1%}")

    report.append("\n## Detection Report by Cell Type\n")
    report.append("| Cell Type | Canonical | Auxiliary | Detection Rate |")
    report.append("|-----------|-----------|-----------|----------------|")

    for cell_type, markers in marker_library.items():
        n_canonical = len(markers["canonical"])
        n_auxiliary = len(markers["auxiliary"])
        detection = detection_report[cell_type]["detection_rate"]

        report.append(
            f"| {cell_type} | {n_canonical} | {n_auxiliary} | {detection:.1%} |"
        )

    if warnings:
        report.append("\n## ⚠️ Warnings\n")
        for warning in warnings:
            report.append(f"### {warning['cell_type']} ({warning['severity']})")
            report.append(f"- Detection rate: {warning['detection_rate']:.1%}")
            report.append(f"- Missing markers: {', '.join(warning['missing_markers'])}")

    report.append("\n## Next Steps")
    report.append("Proceed to Phase 3 (Multi-Strategy Annotation)?")

    return "\n".join(report)
```

## Complete Phase 2 Execution Script

```python
def execute_phase2(adata, annotation_strategy, custom_marker_file=None):
    """
    Main execution function for Phase 2.
    """
    tissue_type = annotation_strategy["tissue_type"]
    species = annotation_strategy["species"]

    print("Phase 2: Marker Gene Library Preparation")
    print("=" * 50)

    # Step 1: Query databases
    print("\n[Step 1/6] Querying marker databases...")
    database_markers = consolidate_database_markers(tissue_type, species)
    print(f"  Retrieved markers for {len(database_markers)} cell types")

    # Step 2: Load custom markers
    print("\n[Step 2/6] Loading custom markers...")
    custom_markers = load_custom_markers(custom_marker_file)
    if custom_markers:
        print(f"  Loaded custom markers for {len(custom_markers)} cell types")

    # Merge sources
    marker_library = merge_marker_sources(database_markers, custom_markers)

    # Step 3: Validate detection
    print("\n[Step 3/6] Validating marker detection...")
    detection_report = calculate_marker_detection(adata, marker_library)
    warnings = flag_critical_missing(detection_report)

    if warnings:
        print(f"  ⚠️  {len(warnings)} cell types have low marker coverage")

    # Step 4: Prioritize markers
    print("\n[Step 4/6] Scoring and ranking markers...")
    specificity_scores = score_marker_specificity(adata, marker_library)
    ranked_markers = rank_markers(marker_library, specificity_scores, detection_report)

    # Step 5: Resolve conflicts
    print("\n[Step 5/6] Resolving marker ambiguities...")
    ambiguous_markers = identify_ambiguous_markers(marker_library)
    combinations = create_marker_combinations(marker_library, ambiguous_markers)
    logic_doc = document_marker_logic(marker_library, ambiguous_markers, combinations)

    print(f"  Found {len(ambiguous_markers)} ambiguous markers")
    print(f"  Created {len(combinations)} marker combinations")

    # Step 6: Export
    print("\n[Step 6/6] Exporting marker library...")
    library = export_marker_library(
        marker_library,
        ranked_markers,
        detection_report,
        ambiguous_markers,
        combinations
    )

    # Save logic document
    with open("docs/03_plan/marker_logic.md", 'w') as f:
        f.write(logic_doc)

    # Generate summary
    summary = generate_summary_report(marker_library, detection_report, warnings)
    print("\n" + summary)

    return library, detection_report, warnings
```

## Quality Gates

Before proceeding to Phase 3, verify:

- [ ] Marker library contains all expected cell types from Phase 1
- [ ] Average detection rate ≥ 70%
- [ ] No cell type has detection rate < 50% (or user acknowledges risk)
- [ ] Ambiguous markers documented with disambiguation logic
- [ ] User approves the marker library

## Output Files

```
docs/03_plan/
├── marker_gene_library.yaml      # Main library file
├── marker_logic.md               # Human-readable logic
├── marker_detection_report.csv   # Detection statistics
└── marker_specificity_scores.csv # Specificity metrics
```

## STOP — Wait for User Approval

<IRON-LAW>
After Phase 2 completes, present the summary report and **STOP**.

Do NOT automatically proceed to Phase 3.
Do NOT begin running annotation algorithms.

The user must review the marker library and explicitly approve before annotation begins.

If warnings exist (low detection rates), the user must acknowledge the risks.
</IRON-LAW>
