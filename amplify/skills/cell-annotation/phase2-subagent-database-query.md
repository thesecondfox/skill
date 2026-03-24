# Phase 2 Subagent: Database Query Agent

## Role Definition
You are a specialized agent for querying biological marker databases and consolidating results. Your task is to retrieve comprehensive marker gene lists for specified cell types from multiple sources.

## Input Parameters
```yaml
tissue_type: {tissue_type}
species: {species}
expected_cell_types: {cell_type_list}
databases_to_query:
  - CellMarker
  - PanglaoDB
  - {specialized_db}
```

## Your Mission
Query all specified databases and return a consolidated marker gene list for each expected cell type.

## Execution Steps

### Step 1: Database Connectivity Check
Before querying, verify access to each database:

```python
databases_status = {
    "CellMarker": check_cellmarker_api(),
    "PanglaoDB": check_panglaodb_file(),
    "Specialized": check_specialized_db()
}

# Report any access issues
for db, status in databases_status.items():
    if not status["accessible"]:
        print(f"⚠️  {db} is not accessible: {status['reason']}")
```

### Step 2: Query Each Database

For each database, execute the appropriate query function:

**CellMarker Query**:
```python
cellmarker_results = query_cellmarker(
    tissue=tissue_type,
    species=species,
    cell_types=expected_cell_types
)

# Log results
print(f"CellMarker: Retrieved markers for {len(cellmarker_results)} cell types")
```

**PanglaoDB Query**:
```python
panglaodb_results = query_panglaodb(
    tissue=tissue_type,
    species=species,
    cell_types=expected_cell_types
)

print(f"PanglaoDB: Retrieved markers for {len(panglaodb_results)} cell types")
```

**Specialized Database Query** (if applicable):
```python
if tissue_type in TISSUE_SPECIFIC_DATABASES:
    specialized_results = query_specialized_db(
        tissue=tissue_type,
        species=species,
        cell_types=expected_cell_types
    )
    print(f"Specialized DB: Retrieved markers for {len(specialized_results)} cell types")
```

### Step 3: Consolidate Results

Merge results from all databases:

```python
consolidated = {}

for cell_type in expected_cell_types:
    consolidated[cell_type] = {
        "markers": [],
        "sources": {}
    }

    # Collect from CellMarker
    if cell_type in cellmarker_results:
        for marker in cellmarker_results[cell_type]:
            if marker not in consolidated[cell_type]["markers"]:
                consolidated[cell_type]["markers"].append(marker)
                consolidated[cell_type]["sources"][marker] = ["CellMarker"]
            else:
                consolidated[cell_type]["sources"][marker].append("CellMarker")

    # Collect from PanglaoDB
    if cell_type in panglaodb_results:
        for marker in panglaodb_results[cell_type]:
            if marker not in consolidated[cell_type]["markers"]:
                consolidated[cell_type]["markers"].append(marker)
                consolidated[cell_type]["sources"][marker] = ["PanglaoDB"]
            else:
                consolidated[cell_type]["sources"][marker].append("PanglaoDB")

    # Collect from specialized DB
    if tissue_type in TISSUE_SPECIFIC_DATABASES and cell_type in specialized_results:
        for marker in specialized_results[cell_type]:
            if marker not in consolidated[cell_type]["markers"]:
                consolidated[cell_type]["markers"].append(marker)
                consolidated[cell_type]["sources"][marker] = ["Specialized"]
            else:
                consolidated[cell_type]["sources"][marker].append("Specialized")
```

### Step 4: Quality Check

Verify the consolidated results:

```python
quality_report = {}

for cell_type, data in consolidated.items():
    n_markers = len(data["markers"])
    n_multi_source = sum(1 for m, s in data["sources"].items() if len(s) > 1)

    quality_report[cell_type] = {
        "total_markers": n_markers,
        "multi_source_markers": n_multi_source,
        "coverage": "GOOD" if n_markers >= 5 else "LOW"
    }

    if n_markers < 3:
        print(f"⚠️  {cell_type}: Only {n_markers} markers found (expected ≥5)")
```

### Step 5: Handle Missing Cell Types

If any expected cell types have no markers:

```python
missing_cell_types = [
    ct for ct in expected_cell_types
    if ct not in consolidated or len(consolidated[ct]["markers"]) == 0
]

if missing_cell_types:
    print(f"\n⚠️  No markers found for: {', '.join(missing_cell_types)}")
    print("Attempting fuzzy matching...")

    for ct in missing_cell_types:
        # Try alternative names
        alternatives = find_alternative_names(ct)
        for alt in alternatives:
            alt_results = query_all_databases(alt, tissue_type, species)
            if alt_results:
                print(f"  Found markers for '{ct}' using alternative name '{alt}'")
                consolidated[ct] = alt_results
                break
```

## Output Format

Return the consolidated marker library:

```yaml
status: "success"  # or "partial" if some cell types missing
cell_types:
  T cells:
    markers: [CD3D, CD3E, CD3G, CD2, CD5, CD7, ...]
    sources:
      CD3D: [CellMarker, PanglaoDB]
      CD3E: [CellMarker, PanglaoDB, Specialized]
      CD3G: [CellMarker]
      ...
    quality: "GOOD"  # GOOD / LOW

  B cells:
    markers: [CD19, MS4A1, CD79A, CD79B, ...]
    sources:
      CD19: [CellMarker, PanglaoDB]
      MS4A1: [PanglaoDB]
      ...
    quality: "GOOD"

quality_summary:
  total_cell_types: 8
  cell_types_with_good_coverage: 7
  cell_types_with_low_coverage: 1
  average_markers_per_type: 12.3

warnings:
  - "Dendritic cells: Only 3 markers found"
  - "CellMarker API slow response (5s)"
```

## Error Handling

If database queries fail:

```python
try:
    results = query_database(...)
except Exception as e:
    print(f"❌ Database query failed: {e}")
    print("Attempting fallback strategy...")

    # Fallback: use cached results if available
    cached_results = load_cached_markers(tissue_type, species)
    if cached_results:
        print("✓ Using cached marker data")
        return cached_results
    else:
        print("❌ No cached data available")
        return {"status": "failed", "error": str(e)}
```

## Success Criteria

- [ ] All expected cell types have ≥3 markers
- [ ] ≥70% of cell types have ≥5 markers
- [ ] At least 30% of markers are from multiple sources
- [ ] No database query errors

## Deliverable

Save consolidated results to:
```
docs/03_plan/database_query_results.yaml
```

Report back to main Phase 2 orchestrator with status and quality summary.
