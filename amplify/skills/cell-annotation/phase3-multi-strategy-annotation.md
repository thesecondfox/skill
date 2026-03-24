# Phase 3: Multi-Strategy Annotation Execution

## Overview
This phase executes the annotation strategy approved in Phase 1, using the marker library prepared in Phase 2. It runs multiple annotation methods in parallel, then integrates results through consensus voting.

## Prerequisites
- Phase 1 completed: `annotation-strategy.yaml` exists
- Phase 2 completed: `marker_gene_library.yaml` exists
- Data loaded with clustering results

## Execution Flow

```
Step 1: Strategy Dispatch
  ├─ Read annotation strategy from Phase 1
  ├─ Identify which methods to run
  └─ Prepare method-specific parameters

Step 2: Parallel Method Execution
  ├─ Method A: Marker-Based Annotation (always)
  ├─ Method B: Reference-Based Annotation (if strategy includes)
  └─ Method C: LLM-Assisted Annotation (for ambiguous clusters)

Step 3: Result Collection
  ├─ Gather annotations from all methods
  ├─ Compute confidence scores
  └─ Flag discrepancies

Step 4: Consensus Integration
  ├─ Cross-method agreement analysis
  ├─ Voting mechanism for final labels
  └─ Confidence scoring

Step 5: Quality Validation
  ├─ Domain sanity check (biological plausibility)
  ├─ Proportion check (vs expected from Phase 1)
  └─ Marker expression validation

Step 6: Manual Review Triage
  ├─ Flag low-confidence clusters
  ├─ Flag conflicting annotations
  └─ Generate review report
```

## Step 1: Strategy Dispatch

### 1.1 Load Strategy Configuration

```python
import yaml
import scanpy as sc

def load_annotation_strategy():
    """
    Load the approved annotation strategy from Phase 1.
    """
    with open("docs/03_plan/annotation-strategy.yaml", 'r') as f:
        strategy = yaml.safe_load(f)

    return strategy

def load_marker_library():
    """
    Load the marker library from Phase 2.
    """
    with open("docs/03_plan/marker_gene_library.yaml", 'r') as f:
        library = yaml.safe_load(f)

    return library

strategy = load_annotation_strategy()
marker_library = load_marker_library()

print(f"Strategy: {strategy['primary_method']}")
print(f"Marker library: {len(marker_library['cell_types'])} cell types")
```

### 1.2 Identify Methods to Execute

```python
def parse_annotation_methods(strategy):
    """
    Extract which annotation methods to run from strategy.
    """
    methods_to_run = []

    if strategy["primary_method"] == "marker-based":
        methods_to_run.append({
            "name": "marker-based",
            "priority": 1,
            "params": strategy["strategy_details"]["step1"]
        })

    elif strategy["primary_method"] == "reference-based":
        methods_to_run.append({
            "name": "reference-based",
            "priority": 1,
            "params": strategy["strategy_details"]["step1"]
        })

    elif strategy["primary_method"] == "hybrid":
        # Hybrid: run both methods
        methods_to_run.append({
            "name": "marker-based",
            "priority": 1,
            "params": strategy["strategy_details"]["step1"]
        })
        methods_to_run.append({
            "name": "reference-based",
            "priority": 2,
            "params": strategy["strategy_details"]["step2"]
        })

    # LLM-assisted for ambiguous clusters (always available)
    if "step3" in strategy["strategy_details"]:
        methods_to_run.append({
            "name": "llm-assisted",
            "priority": 3,
            "params": strategy["strategy_details"]["step3"]
        })

    return methods_to_run

methods = parse_annotation_methods(strategy)
print(f"Methods to execute: {[m['name'] for m in methods]}")
```

## Step 2: Parallel Method Execution

### 2.1 Method A: Marker-Based Annotation

```python
def marker_based_annotation(adata, marker_library, params):
    """
    Annotate clusters using marker gene scoring.

    Algorithm:
    1. For each cluster, compute marker scores for each cell type
    2. Assign cell type with highest score
    3. Compute confidence based on score gap
    """
    import numpy as np
    import pandas as pd

    cluster_key = params.get("cluster_key", "leiden")
    confidence_threshold = params.get("confidence_threshold", 0.7)

    results = {
        "method": "marker-based",
        "annotations": {},
        "confidence_scores": {},
        "marker_scores": {}
    }

    # For each cluster
    for cluster_id in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_cells = adata[cluster_mask]

        # Score each cell type
        cell_type_scores = {}

        for cell_type, markers in marker_library["cell_types"].items():
            canonical = markers["canonical_markers"]
            auxiliary = markers.get("auxiliary_markers", [])

            # Compute score
            score = compute_marker_score(
                cluster_cells,
                canonical,
                auxiliary,
                weights={"canonical": 1.0, "auxiliary": 0.5}
            )

            cell_type_scores[cell_type] = score

        # Assign cell type with highest score
        sorted_scores = sorted(
            cell_type_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )

        best_type, best_score = sorted_scores[0]
        second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else 0

        # Confidence = gap between top 2 scores
        confidence = (best_score - second_score) / (best_score + 1e-6)

        results["annotations"][cluster_id] = best_type
        results["confidence_scores"][cluster_id] = confidence
        results["marker_scores"][cluster_id] = cell_type_scores

        print(f"Cluster {cluster_id}: {best_type} (confidence: {confidence:.2f})")

    return results


def compute_marker_score(cluster_cells, canonical_markers, auxiliary_markers, weights):
    """
    Compute marker score for a cluster.

    Score = (canonical_score * w1 + auxiliary_score * w2) / (w1 + w2)
    """
    import numpy as np

    def score_marker_set(cells, markers):
        """Score a set of markers."""
        if not markers:
            return 0.0

        scores = []
        for marker in markers:
            if marker not in cells.var_names:
                continue

            # Get expression
            expr = cells[:, marker].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()

            # Score = mean expression * % expressing
            mean_expr = np.mean(expr)
            pct_expr = np.sum(expr > 0) / len(expr)

            marker_score = mean_expr * pct_expr
            scores.append(marker_score)

        return np.mean(scores) if scores else 0.0

    canonical_score = score_marker_set(cluster_cells, canonical_markers)
    auxiliary_score = score_marker_set(cluster_cells, auxiliary_markers)

    total_score = (
        canonical_score * weights["canonical"] +
        auxiliary_score * weights["auxiliary"]
    ) / (weights["canonical"] + weights["auxiliary"])

    return total_score
```

### 2.2 Method B: Reference-Based Annotation

```python
def reference_based_annotation(adata, params):
    """
    Annotate using reference dataset (SingleR, CellTypist, etc.).
    """
    tool = params.get("tool", "SingleR")
    reference = params.get("reference", None)
    confidence_threshold = params.get("confidence_threshold", 0.7)

    if tool == "SingleR":
        return run_singler(adata, reference, confidence_threshold)
    elif tool == "CellTypist":
        return run_celltypist(adata, reference, confidence_threshold)
    elif tool == "Azimuth":
        return run_azimuth(adata, reference, confidence_threshold)
    else:
        raise ValueError(f"Unknown reference-based tool: {tool}")


def run_singler(adata, reference_name, confidence_threshold):
    """
    Run SingleR annotation.

    Note: This requires R and SingleR installed.
    We'll use rpy2 to call R from Python.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    # Load SingleR in R
    ro.r('''
        library(SingleR)
        library(celldex)
    ''')

    # Load reference
    if reference_name == "Monaco_Immune_2019":
        ro.r('ref <- MonacoImmuneData()')
    elif reference_name == "HPCA_2018":
        ro.r('ref <- HumanPrimaryCellAtlasData()')
    else:
        raise ValueError(f"Unknown reference: {reference_name}")

    # Convert adata to R format
    # (Simplified - in practice, use anndata2ri)
    expr_matrix = adata.X.T  # Genes × Cells

    ro.globalenv['test_data'] = expr_matrix
    ro.globalenv['gene_names'] = list(adata.var_names)

    # Run SingleR
    ro.r('''
        rownames(test_data) <- gene_names
        pred <- SingleR(
            test = test_data,
            ref = ref,
            labels = ref$label.main,
            de.method = "wilcox"
        )
    ''')

    # Get results
    predictions = ro.r('pred$labels')
    scores = ro.r('pred$scores')

    # Convert to cluster-level annotations
    cluster_key = "leiden"
    results = {
        "method": "reference-based",
        "tool": "SingleR",
        "reference": reference_name,
        "annotations": {},
        "confidence_scores": {}
    }

    for cluster_id in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_preds = predictions[cluster_mask]

        # Majority vote
        from collections import Counter
        vote_counts = Counter(cluster_preds)
        best_type, count = vote_counts.most_common(1)[0]

        # Confidence = proportion of cells agreeing
        confidence = count / len(cluster_preds)

        results["annotations"][cluster_id] = best_type
        results["confidence_scores"][cluster_id] = confidence

    return results


def run_celltypist(adata, model_name, confidence_threshold):
    """
    Run CellTypist annotation.

    CellTypist is a Python package - easier to use than SingleR.
    """
    import celltypist
    from celltypist import models

    # Download model if needed
    models.download_models(model=model_name)

    # Load model
    model = models.Model.load(model=model_name)

    # Run prediction
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True  # Cluster-level voting
    )

    # Extract results
    cluster_key = "leiden"
    results = {
        "method": "reference-based",
        "tool": "CellTypist",
        "reference": model_name,
        "annotations": {},
        "confidence_scores": {}
    }

    for cluster_id in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster_id

        # Get majority-voted label
        label = predictions.predicted_labels.loc[cluster_mask, "majority_voting"].iloc[0]

        # Get confidence (proportion of cells with this label)
        cluster_labels = predictions.predicted_labels.loc[cluster_mask, "predicted_labels"]
        confidence = (cluster_labels == label).sum() / len(cluster_labels)

        results["annotations"][cluster_id] = label
        results["confidence_scores"][cluster_id] = confidence

    return results
```

### 2.3 Method C: LLM-Assisted Annotation

```python
def llm_assisted_annotation(adata, marker_library, ambiguous_clusters, params):
    """
    Use LLM (GPT-4/Claude) to interpret marker genes for ambiguous clusters.
    """
    import anthropic

    client = anthropic.Anthropic()

    results = {
        "method": "llm-assisted",
        "annotations": {},
        "confidence_scores": {},
        "reasoning": {}
    }

    for cluster_id in ambiguous_clusters:
        # Get top marker genes for this cluster
        cluster_mask = adata.obs["leiden"] == cluster_id
        cluster_cells = adata[cluster_mask]

        # Compute differential expression
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster_id],
            reference="rest",
            method="wilcoxon"
        )

        top_markers = adata.uns["rank_genes_groups"]["names"][cluster_id][:20]

        # Prepare prompt for LLM
        prompt = f"""
You are a single-cell biologist. Based on the following top 20 marker genes
for a cluster in {strategy['tissue_type']} tissue, identify the most likely
cell type.

Top marker genes (ranked by differential expression):
{', '.join(top_markers)}

Expected cell types in this tissue:
{', '.join(marker_library['cell_types'].keys())}

Provide:
1. Most likely cell type
2. Confidence level (HIGH/MEDIUM/LOW)
3. Brief reasoning (2-3 sentences)

Format your response as:
Cell Type: [type]
Confidence: [HIGH/MEDIUM/LOW]
Reasoning: [explanation]
"""

        # Call LLM
        message = client.messages.create(
            model="claude-opus-4-6",
            max_tokens=500,
            messages=[{"role": "user", "content": prompt}]
        )

        response_text = message.content[0].text

        # Parse response
        cell_type = extract_field(response_text, "Cell Type")
        confidence_str = extract_field(response_text, "Confidence")
        reasoning = extract_field(response_text, "Reasoning")

        # Convert confidence to numeric
        confidence_map = {"HIGH": 0.9, "MEDIUM": 0.7, "LOW": 0.5}
        confidence = confidence_map.get(confidence_str, 0.5)

        results["annotations"][cluster_id] = cell_type
        results["confidence_scores"][cluster_id] = confidence
        results["reasoning"][cluster_id] = reasoning

        print(f"Cluster {cluster_id}: {cell_type} ({confidence_str})")
        print(f"  Reasoning: {reasoning}")

    return results


def extract_field(text, field_name):
    """Extract field value from LLM response."""
    import re
    pattern = f"{field_name}:\\s*(.+?)(?:\\n|$)"
    match = re.search(pattern, text)
    return match.group(1).strip() if match else "Unknown"
```

### 2.4 Parallel Dispatch

```python
def execute_annotation_methods_parallel(adata, marker_library, strategy, methods):
    """
    Dispatch all annotation methods in parallel using Task tool.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results = {}

    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {}

        for method_config in methods:
            method_name = method_config["name"]
            params = method_config["params"]

            if method_name == "marker-based":
                future = executor.submit(
                    marker_based_annotation,
                    adata, marker_library, params
                )
                futures[future] = method_name

            elif method_name == "reference-based":
                future = executor.submit(
                    reference_based_annotation,
                    adata, params
                )
                futures[future] = method_name

            elif method_name == "llm-assisted":
                # Only run on ambiguous clusters
                ambiguous = identify_ambiguous_clusters(adata, results)
                future = executor.submit(
                    llm_assisted_annotation,
                    adata, marker_library, ambiguous, params
                )
                futures[future] = method_name

        # Collect results
        for future in as_completed(futures):
            method_name = futures[future]
            try:
                result = future.result()
                results[method_name] = result
                print(f"✓ {method_name} completed")
            except Exception as e:
                print(f"✗ {method_name} failed: {e}")
                results[method_name] = {"error": str(e)}

    return results
```

## Step 3: Result Collection

```python
def collect_annotation_results(results):
    """
    Consolidate results from all methods into a unified format.
    """
    consolidated = {
        "methods_run": list(results.keys()),
        "clusters": {},
        "method_agreement": {}
    }

    # Get all cluster IDs
    all_clusters = set()
    for method_name, method_results in results.items():
        if "annotations" in method_results:
            all_clusters.update(method_results["annotations"].keys())

    # For each cluster, collect annotations from all methods
    for cluster_id in all_clusters:
        consolidated["clusters"][cluster_id] = {
            "annotations": {},
            "confidence_scores": {},
            "consensus": None,
            "consensus_confidence": None
        }

        for method_name, method_results in results.items():
            if "annotations" in method_results:
                annotation = method_results["annotations"].get(cluster_id)
                confidence = method_results["confidence_scores"].get(cluster_id)

                if annotation:
                    consolidated["clusters"][cluster_id]["annotations"][method_name] = annotation
                    consolidated["clusters"][cluster_id]["confidence_scores"][method_name] = confidence

    return consolidated
```

## Step 4: Consensus Integration

```python
def compute_consensus_annotations(consolidated, strategy):
    """
    Integrate annotations from multiple methods using voting.
    """
    validation_config = strategy["strategy_details"]["validation"]
    agreement_threshold = validation_config["cross_method_agreement_threshold"]

    for cluster_id, cluster_data in consolidated["clusters"].items():
        annotations = cluster_data["annotations"]
        confidences = cluster_data["confidence_scores"]

        if not annotations:
            cluster_data["consensus"] = "Unknown"
            cluster_data["consensus_confidence"] = 0.0
            cluster_data["agreement_status"] = "NO_ANNOTATIONS"
            continue

        # Count votes (weighted by confidence)
        from collections import defaultdict
        weighted_votes = defaultdict(float)

        for method, annotation in annotations.items():
            confidence = confidences[method]
            weighted_votes[annotation] += confidence

        # Get top annotation
        sorted_votes = sorted(
            weighted_votes.items(),
            key=lambda x: x[1],
            reverse=True
        )

        consensus_type, consensus_weight = sorted_votes[0]
        total_weight = sum(weighted_votes.values())

        # Consensus confidence = proportion of weighted votes
        consensus_confidence = consensus_weight / total_weight

        # Check agreement
        if len(annotations) > 1:
            # All methods agree?
            if len(set(annotations.values())) == 1:
                agreement_status = "FULL_AGREEMENT"
            elif consensus_confidence >= agreement_threshold:
                agreement_status = "MAJORITY_AGREEMENT"
            else:
                agreement_status = "DISAGREEMENT"
        else:
            agreement_status = "SINGLE_METHOD"

        cluster_data["consensus"] = consensus_type
        cluster_data["consensus_confidence"] = consensus_confidence
        cluster_data["agreement_status"] = agreement_status

        print(f"Cluster {cluster_id}: {consensus_type} "
              f"(confidence: {consensus_confidence:.2f}, {agreement_status})")

    return consolidated
```

## Step 5: Quality Validation

### 5.1 Domain Sanity Check

```python
def domain_sanity_check(consolidated, marker_library, strategy):
    """
    Check if annotations are biologically plausible.
    """
    tissue_type = strategy["tissue_type"]
    expected_cell_types = strategy["expected_cell_types"]

    sanity_report = {
        "passed": [],
        "warnings": [],
        "failures": []
    }

    # Check 1: Are all annotations in expected cell types?
    for cluster_id, cluster_data in consolidated["clusters"].items():
        consensus = cluster_data["consensus"]

        if consensus not in [ct["name"] for ct in expected_cell_types]:
            sanity_report["warnings"].append({
                "cluster": cluster_id,
                "issue": "UNEXPECTED_CELL_TYPE",
                "annotation": consensus,
                "message": f"Cell type '{consensus}' not in expected types for {tissue_type}"
            })

    # Check 2: Cell type proportions
    total_cells = sum(
        len(consolidated["clusters"][cid]["annotations"])
        for cid in consolidated["clusters"]
    )

    for expected_ct in expected_cell_types:
        ct_name = expected_ct["name"]
        expected_range = expected_ct["expected_proportion"]

        # Count cells of this type
        ct_clusters = [
            cid for cid, data in consolidated["clusters"].items()
            if data["consensus"] == ct_name
        ]

        # Calculate proportion (simplified - should count actual cells)
        ct_proportion = len(ct_clusters) / len(consolidated["clusters"])

        # Parse expected range (e.g., "60-70%")
        import re
        match = re.match(r"(\d+)-(\d+)%", expected_range)
        if match:
            min_pct, max_pct = int(match.group(1)), int(match.group(2))

            if not (min_pct/100 <= ct_proportion <= max_pct/100):
                sanity_report["warnings"].append({
                    "cell_type": ct_name,
                    "issue": "PROPORTION_MISMATCH",
                    "expected": expected_range,
                    "observed": f"{ct_proportion*100:.1f}%",
                    "message": f"{ct_name} proportion outside expected range"
                })

    # Check 3: Marker expression validation
    # For each annotation, verify key markers are actually expressed
    for cluster_id, cluster_data in consolidated["clusters"].items():
        consensus = cluster_data["consensus"]

        if consensus in marker_library["cell_types"]:
            canonical_markers = marker_library["cell_types"][consensus]["canonical_markers"]

            # Check if canonical markers are expressed in this cluster
            # (Implementation depends on having cluster expression data)
            # This is a placeholder
            pass

    return sanity_report
```

### 5.2 Confidence-Based Validation

```python
def validate_confidence_levels(consolidated, strategy):
    """
    Flag clusters with low confidence or disagreement.
    """
    validation_config = strategy["strategy_details"]["validation"]
    manual_review_triggers = validation_config["manual_review_triggers"]

    review_list = []

    for cluster_id, cluster_data in consolidated["clusters"].items():
        consensus_confidence = cluster_data["consensus_confidence"]
        agreement_status = cluster_data["agreement_status"]

        # Check triggers
        needs_review = False
        reasons = []

        # Trigger 1: Low confidence
        if consensus_confidence < 0.7:
            needs_review = True
            reasons.append(f"Low confidence ({consensus_confidence:.2f})")

        # Trigger 2: Methods disagree
        if agreement_status == "DISAGREEMENT":
            needs_review = True
            reasons.append("Methods disagree")

        # Trigger 3: Small cluster size
        # (Would need actual cell counts)

        if needs_review:
            review_list.append({
                "cluster_id": cluster_id,
                "consensus": cluster_data["consensus"],
                "confidence": consensus_confidence,
                "agreement": agreement_status,
                "reasons": reasons,
                "annotations_by_method": cluster_data["annotations"]
            })

    return review_list
```

## Step 6: Manual Review Triage

```python
def generate_review_report(consolidated, review_list, sanity_report):
    """
    Generate a report for manual review of flagged clusters.
    """
    report = []

    report.append("# Manual Review Required\n")
    report.append(f"Total clusters: {len(consolidated['clusters'])}")
    report.append(f"Flagged for review: {len(review_list)}\n")

    report.append("## High-Priority Review\n")

    for item in review_list:
        report.append(f"### Cluster {item['cluster_id']}")
        report.append(f"**Consensus**: {item['consensus']}")
        report.append(f"**Confidence**: {item['confidence']:.2f}")
        report.append(f"**Agreement**: {item['agreement']}")
        report.append(f"**Reasons**: {', '.join(item['reasons'])}\n")

        report.append("**Annotations by method**:")
        for method, annotation in item['annotations_by_method'].items():
            report.append(f"- {method}: {annotation}")
        report.append("")

    report.append("\n## Sanity Check Warnings\n")

    for warning in sanity_report["warnings"]:
        report.append(f"- **{warning['issue']}**: {warning['message']}")

    return "\n".join(report)
```

## Complete Phase 3 Execution

```python
def execute_phase3(adata, strategy_path, marker_library_path):
    """
    Main execution function for Phase 3.
    """
    print("Phase 3: Multi-Strategy Annotation Execution")
    print("=" * 50)

    # Step 1: Load configuration
    print("\n[Step 1/6] Loading strategy and marker library...")
    strategy = load_annotation_strategy()
    marker_library = load_marker_library()

    # Step 2: Parse methods
    print("\n[Step 2/6] Parsing annotation methods...")
    methods = parse_annotation_methods(strategy)
    print(f"  Methods to run: {[m['name'] for m in methods]}")

    # Step 3: Execute methods in parallel
    print("\n[Step 3/6] Executing annotation methods...")
    results = execute_annotation_methods_parallel(
        adata, marker_library, strategy, methods
    )

    # Step 4: Collect and consolidate
    print("\n[Step 4/6] Consolidating results...")
    consolidated = collect_annotation_results(results)
    consolidated = compute_consensus_annotations(consolidated, strategy)

    # Step 5: Quality validation
    print("\n[Step 5/6] Running quality validation...")
    sanity_report = domain_sanity_check(consolidated, marker_library, strategy)
    review_list = validate_confidence_levels(consolidated, strategy)

    print(f"  Sanity check: {len(sanity_report['warnings'])} warnings")
    print(f"  Manual review: {len(review_list)} clusters flagged")

    # Step 6: Generate reports
    print("\n[Step 6/6] Generating reports...")
    review_report = generate_review_report(consolidated, review_list, sanity_report)

    # Save results
    save_annotation_results(consolidated, review_list, sanity_report)

    # Print summary
    print_phase3_summary(consolidated, review_list, sanity_report)

    return consolidated, review_list, sanity_report
```

## Output Files

```
docs/04_annotation/
├── annotation_results.yaml          # All annotations with confidence
├── consensus_annotations.csv        # Final consensus labels
├── method_comparison.csv            # Cross-method comparison
├── manual_review_list.yaml          # Clusters needing review
├── sanity_check_report.md           # Biological plausibility check
└── annotation_summary.md            # Human-readable summary
```

## Quality Gates

Before proceeding to Phase 4, verify:

- [ ] All clusters have consensus annotations
- [ ] ≥80% of clusters have confidence >0.7
- [ ] No critical sanity check failures
- [ ] Manual review list is reasonable (<20% of clusters)
- [ ] User reviews and approves flagged clusters

## STOP — Wait for User Approval

<IRON-LAW>
After Phase 3 completes, present the annotation summary and **STOP**.

Do NOT automatically proceed to Phase 4.
Do NOT apply annotations to the data object yet.

The user must review:
1. Consensus annotations
2. Flagged clusters requiring manual review
3. Sanity check warnings

User must explicitly approve before finalizing annotations.
</IRON-LAW>
