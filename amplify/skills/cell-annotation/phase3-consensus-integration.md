# Phase 3 Orchestration: Consensus Integration

## Overview
This document describes how to integrate results from multiple annotation methods (marker-based, reference-based, LLM-assisted) into a final consensus annotation with confidence scores.

## Input
Results from all annotation methods:
```yaml
marker_based_results: {annotations, confidence_scores, evidence}
reference_based_results: {annotations, confidence_scores, validation}
llm_assisted_results: {annotations, confidence_scores, reasoning}  # optional
```

## Consensus Algorithm

### Step 1: Collect All Annotations

```python
def collect_all_annotations(results_dict):
    """
    Consolidate annotations from all methods.
    """
    all_clusters = set()

    # Collect all cluster IDs
    for method_name, results in results_dict.items():
        if "annotations" in results:
            all_clusters.update(results["annotations"].keys())

    # Build consolidated structure
    consolidated = {}

    for cluster_id in all_clusters:
        consolidated[cluster_id] = {
            "methods": {},
            "consensus": None,
            "consensus_confidence": None,
            "agreement_status": None
        }

        for method_name, results in results_dict.items():
            if cluster_id in results.get("annotations", {}):
                consolidated[cluster_id]["methods"][method_name] = {
                    "annotation": results["annotations"][cluster_id],
                    "confidence": results["confidence_scores"][cluster_id]
                }

    return consolidated
```

### Step 2: Voting Mechanism

```python
def compute_consensus_vote(cluster_data, voting_strategy="weighted"):
    """
    Compute consensus annotation using voting.

    Strategies:
    - weighted: votes weighted by confidence
    - majority: simple majority vote
    - hierarchical: marker-based > reference-based > LLM
    """
    methods = cluster_data["methods"]

    if not methods:
        return None, 0.0, "NO_ANNOTATIONS"

    if voting_strategy == "weighted":
        return weighted_voting(methods)
    elif voting_strategy == "majority":
        return majority_voting(methods)
    elif voting_strategy == "hierarchical":
        return hierarchical_voting(methods)
    else:
        raise ValueError(f"Unknown voting strategy: {voting_strategy}")


def weighted_voting(methods):
    """
    Weighted voting: each method's vote weighted by its confidence.
    """
    from collections import defaultdict

    weighted_votes = defaultdict(float)

    for method_name, data in methods.items():
        annotation = data["annotation"]
        confidence = data["confidence"]

        # Weight by confidence
        weighted_votes[annotation] += confidence

    # Get winner
    sorted_votes = sorted(
        weighted_votes.items(),
        key=lambda x: x[1],
        reverse=True
    )

    if not sorted_votes:
        return None, 0.0, "NO_VOTES"

    winner, winner_weight = sorted_votes[0]
    total_weight = sum(weighted_votes.values())

    # Consensus confidence = proportion of weighted votes
    consensus_confidence = winner_weight / total_weight

    # Agreement status
    if len(set(m["annotation"] for m in methods.values())) == 1:
        agreement = "FULL_AGREEMENT"
    elif consensus_confidence >= 0.8:
        agreement = "STRONG_MAJORITY"
    elif consensus_confidence >= 0.6:
        agreement = "WEAK_MAJORITY"
    else:
        agreement = "DISAGREEMENT"

    return winner, consensus_confidence, agreement


def majority_voting(methods):
    """
    Simple majority vote (unweighted).
    """
    from collections import Counter

    votes = [m["annotation"] for m in methods.values()]
    vote_counts = Counter(votes)

    winner, count = vote_counts.most_common(1)[0]
    consensus_confidence = count / len(votes)

    if len(vote_counts) == 1:
        agreement = "FULL_AGREEMENT"
    elif count > len(votes) / 2:
        agreement = "MAJORITY"
    else:
        agreement = "DISAGREEMENT"

    return winner, consensus_confidence, agreement


def hierarchical_voting(methods):
    """
    Hierarchical: prefer marker-based > reference-based > LLM.
    """
    priority_order = ["marker-based", "reference-based", "llm-assisted"]

    for method_name in priority_order:
        if method_name in methods:
            data = methods[method_name]
            return data["annotation"], data["confidence"], f"HIERARCHICAL_{method_name.upper()}"

    return None, 0.0, "NO_METHODS"
```

### Step 3: Conflict Resolution

```python
def resolve_conflicts(consolidated, marker_library, adata):
    """
    For clusters with disagreement, apply tie-breaking rules.
    """
    conflicts = []

    for cluster_id, data in consolidated.items():
        if data["agreement_status"] == "DISAGREEMENT":
            conflicts.append(cluster_id)

            # Tie-breaking strategy
            methods = data["methods"]

            # Rule 1: If marker-based has high confidence, trust it
            if "marker-based" in methods:
                mb_conf = methods["marker-based"]["confidence"]
                if mb_conf >= 0.8:
                    data["consensus"] = methods["marker-based"]["annotation"]
                    data["consensus_confidence"] = mb_conf
                    data["agreement_status"] = "RESOLVED_BY_MARKER_CONFIDENCE"
                    continue

            # Rule 2: Check marker expression for each candidate
            candidates = {m["annotation"] for m in methods.values()}

            marker_scores = {}
            for candidate in candidates:
                if candidate in marker_library["cell_types"]:
                    markers = marker_library["cell_types"][candidate]["canonical_markers"]
                    score = compute_marker_expression_score(
                        adata, cluster_id, markers
                    )
                    marker_scores[candidate] = score

            # Choose candidate with highest marker score
            if marker_scores:
                best_candidate = max(marker_scores, key=marker_scores.get)
                data["consensus"] = best_candidate
                data["consensus_confidence"] = marker_scores[best_candidate]
                data["agreement_status"] = "RESOLVED_BY_MARKER_EXPRESSION"
            else:
                # Rule 3: Default to weighted vote
                data["consensus"] = data["methods"]["marker-based"]["annotation"]
                data["consensus_confidence"] = 0.5
                data["agreement_status"] = "UNRESOLVED_CONFLICT"

    return conflicts


def compute_marker_expression_score(adata, cluster_id, markers):
    """
    Compute marker expression score for a cluster.
    """
    cluster_mask = adata.obs["leiden"] == cluster_id
    cluster_cells = adata[cluster_mask]

    scores = []
    for marker in markers:
        if marker in adata.var_names:
            expr = cluster_cells[:, marker].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()

            mean_expr = np.mean(expr)
            pct_expr = np.sum(expr > 0) / len(expr)

            scores.append(mean_expr * pct_expr)

    return np.mean(scores) if scores else 0.0
```

### Step 4: Confidence Calibration

```python
def calibrate_confidence(consolidated):
    """
    Adjust confidence scores based on agreement and evidence strength.
    """
    for cluster_id, data in consolidated.items():
        base_confidence = data["consensus_confidence"]
        agreement = data["agreement_status"]

        # Calibration factors
        if agreement == "FULL_AGREEMENT":
            calibrated = min(base_confidence * 1.1, 1.0)  # Boost by 10%
        elif agreement in ["STRONG_MAJORITY", "RESOLVED_BY_MARKER_CONFIDENCE"]:
            calibrated = base_confidence  # No change
        elif agreement in ["WEAK_MAJORITY", "RESOLVED_BY_MARKER_EXPRESSION"]:
            calibrated = base_confidence * 0.9  # Reduce by 10%
        elif agreement in ["DISAGREEMENT", "UNRESOLVED_CONFLICT"]:
            calibrated = base_confidence * 0.7  # Reduce by 30%
        else:
            calibrated = base_confidence

        data["calibrated_confidence"] = calibrated

    return consolidated
```

### Step 5: Generate Consensus Report

```python
def generate_consensus_report(consolidated):
    """
    Create human-readable consensus report.
    """
    report = []

    report.append("# Consensus Annotation Report\n")

    # Summary statistics
    total = len(consolidated)
    full_agreement = sum(
        1 for d in consolidated.values()
        if d["agreement_status"] == "FULL_AGREEMENT"
    )
    high_conf = sum(
        1 for d in consolidated.values()
        if d["calibrated_confidence"] >= 0.8
    )
    low_conf = sum(
        1 for d in consolidated.values()
        if d["calibrated_confidence"] < 0.6
    )

    report.append(f"## Summary")
    report.append(f"- Total clusters: {total}")
    report.append(f"- Full agreement: {full_agreement} ({full_agreement/total*100:.1f}%)")
    report.append(f"- High confidence (≥0.8): {high_conf} ({high_conf/total*100:.1f}%)")
    report.append(f"- Low confidence (<0.6): {low_conf} ({low_conf/total*100:.1f}%)\n")

    # Per-cluster details
    report.append("## Cluster Annotations\n")
    report.append("| Cluster | Consensus | Confidence | Agreement | Methods |")
    report.append("|---------|-----------|------------|-----------|---------|")

    for cluster_id in sorted(consolidated.keys()):
        data = consolidated[cluster_id]
        consensus = data["consensus"]
        conf = data["calibrated_confidence"]
        agreement = data["agreement_status"]

        methods_str = ", ".join(data["methods"].keys())

        report.append(
            f"| {cluster_id} | {consensus} | {conf:.2f} | {agreement} | {methods_str} |"
        )

    # Conflicts and warnings
    conflicts = [
        cid for cid, d in consolidated.items()
        if "CONFLICT" in d["agreement_status"] or "DISAGREEMENT" in d["agreement_status"]
    ]

    if conflicts:
        report.append("\n## ⚠️ Conflicts Requiring Review\n")
        for cluster_id in conflicts:
            data = consolidated[cluster_id]
            report.append(f"### Cluster {cluster_id}")
            report.append(f"**Consensus**: {data['consensus']} (confidence: {data['calibrated_confidence']:.2f})")
            report.append(f"**Status**: {data['agreement_status']}\n")

            report.append("**Annotations by method**:")
            for method, mdata in data["methods"].items():
                report.append(f"- {method}: {mdata['annotation']} (conf: {mdata['confidence']:.2f})")
            report.append("")

    return "\n".join(report)
```

## Complete Consensus Integration

```python
def integrate_consensus(results_dict, marker_library, adata, voting_strategy="weighted"):
    """
    Main function to integrate all annotation results into consensus.
    """
    print("Integrating consensus annotations...")

    # Step 1: Collect
    consolidated = collect_all_annotations(results_dict)
    print(f"  Collected annotations for {len(consolidated)} clusters")

    # Step 2: Vote
    for cluster_id, data in consolidated.items():
        consensus, confidence, agreement = compute_consensus_vote(
            data, voting_strategy
        )
        data["consensus"] = consensus
        data["consensus_confidence"] = confidence
        data["agreement_status"] = agreement

    # Step 3: Resolve conflicts
    conflicts = resolve_conflicts(consolidated, marker_library, adata)
    print(f"  Resolved {len(conflicts)} conflicts")

    # Step 4: Calibrate confidence
    consolidated = calibrate_confidence(consolidated)

    # Step 5: Generate report
    report = generate_consensus_report(consolidated)

    return consolidated, report
```

## Output Format

```yaml
consensus_annotations:
  "0":
    consensus: "CD8+ T cells"
    consensus_confidence: 0.85
    calibrated_confidence: 0.87
    agreement_status: "FULL_AGREEMENT"
    methods:
      marker-based:
        annotation: "CD8+ T cells"
        confidence: 0.85
      reference-based:
        annotation: "CD8+ T cells"
        confidence: 0.92

  "1":
    consensus: "B cells"
    consensus_confidence: 0.78
    calibrated_confidence: 0.78
    agreement_status: "STRONG_MAJORITY"
    methods:
      marker-based:
        annotation: "B cells"
        confidence: 0.72
      reference-based:
        annotation: "B cells"
        confidence: 0.85

  "5":
    consensus: "Dendritic cells"
    consensus_confidence: 0.55
    calibrated_confidence: 0.39
    agreement_status: "DISAGREEMENT"
    methods:
      marker-based:
        annotation: "Dendritic cells"
        confidence: 0.45
      reference-based:
        annotation: "Monocytes"
        confidence: 0.65
      llm-assisted:
        annotation: "Dendritic cells"
        confidence: 0.70

summary:
  total_clusters: 15
  full_agreement: 10
  high_confidence: 11
  low_confidence: 2
  conflicts_resolved: 3
  conflicts_unresolved: 1
```

## Deliverables

Save to:
```
docs/04_annotation/
├── consensus_annotations.yaml
├── consensus_report.md
└── conflict_resolution_log.csv
```
