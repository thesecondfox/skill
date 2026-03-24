# Phase 1 Multi-Agent Orchestration — Main Dispatcher

## Overview
This document describes how to orchestrate the three agents in Phase 1 (Annotation Strategy Selection) using the multi-round deliberation protocol.

## Agent Composition

| Agent | Role | Optimization Target | Prompt Template |
|-------|------|-------------------|-----------------|
| **Agent 1** | Tissue Biology Expert | Biological plausibility and expected cell types | `phase1-agent1-tissue-expert.md` |
| **Agent 2** | Annotation Method Specialist | Optimal method selection and execution plan | `phase1-agent2-method-advisor.md` |
| **Agent 3** | Quality Controller | Data quality and risk mitigation | `phase1-agent3-quality-controller.md` |

## Execution Protocol

### Pre-Dispatch: Context Package Preparation

Before dispatching agents, prepare the shared context:

```python
context_package = {
    # Data characteristics
    "tissue_type": adata.obs['tissue'][0] if 'tissue' in adata.obs else "unknown",
    "species": "human",  # or extract from data
    "n_cells": adata.n_obs,
    "n_genes": adata.n_vars,
    "n_clusters": len(adata.obs['leiden'].unique()),
    "cluster_sizes": adata.obs['leiden'].value_counts().to_dict(),

    # Quality metrics
    "qc_metrics_summary": {
        "median_genes_per_cell": np.median(adata.obs['n_genes']),
        "median_umi_per_cell": np.median(adata.obs['n_counts']),
        "pct_mt_median": np.median(adata.obs['pct_counts_mt']),
    },

    # Clustering quality
    "silhouette_score": compute_silhouette(adata),
    "cluster_purity": compute_purity(adata),

    # Preliminary markers (top 5 per cluster)
    "top_markers_preview": get_top_markers(adata, n_genes=5),

    # Available resources
    "available_references": list_available_references(tissue_type),
    "available_marker_dbs": ["CellMarker", "PanglaoDB", "custom"],
    "compute_resources": "standard",  # or "high-memory", "GPU"

    # User context
    "user_context": user_provided_context,
    "user_preferred_methods": user_preferences.get("methods", None),
    "time_budget": user_preferences.get("time_budget", "standard"),
}
```

### Round 1: Parallel Dispatch

Dispatch all three agents simultaneously with the same context package:

```python
# Agent 1: Tissue Biology Expert
task_1 = Task(
    description="Tissue Biology Expert — evaluate biological context",
    prompt=fill_template(
        "phase1-agent1-tissue-expert.md",
        context_package
    ),
    subagent_type="general-purpose"
)

# Agent 2: Annotation Method Specialist
task_2 = Task(
    description="Method Specialist — recommend annotation strategy",
    prompt=fill_template(
        "phase1-agent2-method-advisor.md",
        context_package
    ),
    subagent_type="general-purpose"
)

# Agent 3: Quality Controller
task_3 = Task(
    description="Quality Controller — audit data and strategy",
    prompt=fill_template(
        "phase1-agent3-quality-controller.md",
        context_package
    ),
    subagent_type="general-purpose"
)

# Wait for all agents to complete
results_round1 = await_all([task_1, task_2, task_3])
```

### Round 1 Synthesis

After Round 1, synthesize the feedback:

```python
synthesis_round1 = {
    "agent1_verdict": extract_verdict(results_round1[0]),
    "agent2_verdict": extract_verdict(results_round1[1]),
    "agent3_verdict": extract_verdict(results_round1[2]),

    "consensus_points": find_agreements(results_round1),
    "conflicts": find_disagreements(results_round1),

    "action_items": {
        "literature_searches": extract_search_queries(results_round1),
        "data_fixes": extract_required_fixes(results_round1),
        "strategy_refinements": extract_recommendations(results_round1),
    }
}
```

### Convergence Check

After each round, check if agents have converged:

```python
def check_convergence(results):
    """
    Convergence criteria:
    1. All agents return PASS verdict
    2. No conflicting recommendations
    3. No new literature searches requested
    4. No blocking issues identified
    """
    verdicts = [extract_verdict(r) for r in results]

    if all(v == "PASS" for v in verdicts):
        conflicts = find_disagreements(results)
        if len(conflicts) == 0:
            return True, "All agents agree — strategy approved"

    if all(v in ["PASS", "CONDITIONAL"] for v in verdicts):
        # Check if conditions are resolvable
        conditions = extract_conditions(results)
        if all(is_resolvable(c) for c in conditions):
            return True, "Conditional approval — address minor issues"

    return False, "Agents have not converged — continue deliberation"
```

### Round 2+: Iterative Refinement

If not converged, address issues and re-dispatch:

```python
# Execute action items from Round 1
if synthesis_round1["action_items"]["literature_searches"]:
    new_literature = execute_searches(
        synthesis_round1["action_items"]["literature_searches"]
    )
    context_package["new_literature"] = new_literature

if synthesis_round1["action_items"]["data_fixes"]:
    apply_fixes(synthesis_round1["action_items"]["data_fixes"])
    # Recompute QC metrics
    context_package["qc_metrics_summary"] = recompute_qc(adata)

# Update context with Round 1 feedback
context_package["round1_feedback"] = {
    "agent1_summary": results_round1[0]["summary"],
    "agent2_summary": results_round1[1]["summary"],
    "agent3_summary": results_round1[2]["summary"],
    "synthesis": synthesis_round1,
}

# Re-dispatch all agents with updated context
results_round2 = dispatch_all_agents(context_package)
```

### Maximum Rounds: 5

Follow the `multi-round-deliberation` protocol:
- Max 5 rounds
- All agents participate in every round
- Each round sees full history of previous rounds
- Stop early if converged

### Escalation Path

If agents don't converge after 5 rounds:

```python
if round_number >= 5 and not converged:
    # Present options to user
    present_to_user({
        "status": "Agents did not reach consensus after 5 rounds",
        "key_disagreements": extract_conflicts(all_results),
        "options": [
            {
                "option": "A",
                "description": "Proceed with majority recommendation",
                "recommendation": get_majority_view(all_results),
                "risks": list_risks(get_majority_view(all_results)),
            },
            {
                "option": "B",
                "description": "Adopt conservative approach",
                "recommendation": get_most_conservative_view(all_results),
                "risks": list_risks(get_most_conservative_view(all_results)),
            },
            {
                "option": "C",
                "description": "Return to Phase 0 (data preprocessing)",
                "reason": "Fundamental data quality issues detected",
            },
            {
                "option": "D",
                "description": "Manual strategy specification",
                "prompt": "Please specify your preferred annotation approach",
            },
        ]
    })
```

## Output Format

After convergence (or user decision), produce the final Phase 1 deliverable:

```yaml
# annotation-strategy.yaml

phase1_status: "converged"  # or "user_override"
rounds_completed: 3
convergence_reason: "All agents approved hybrid strategy"

# Biological context
tissue_type: "PBMC"
species: "human"
expected_cell_types:
  - name: "T cells"
    expected_proportion: "60-70%"
    key_markers: ["CD3D", "CD3E", "CD3G"]
    rarity: "common"
  - name: "B cells"
    expected_proportion: "10-15%"
    key_markers: ["CD19", "MS4A1", "CD79A"]
    rarity: "common"
  # ... more cell types

# Quality assessment
data_quality:
  overall: "good"
  issues_resolved:
    - "Removed 2 low-quality clusters (clusters 8, 12)"
    - "Applied doublet filtering (removed 3.2% of cells)"
  remaining_concerns:
    - "Cluster 5 has moderate doublet score (0.35) — flag for review"

clustering_quality:
  overall: "optimal"
  n_clusters: 15
  resolution: 0.8
  silhouette_score: 0.42

# Annotation strategy (approved by all agents)
primary_method: "hybrid"
strategy_details:
  step1:
    method: "marker-based"
    tool: "scanpy.tl.rank_genes_groups"
    target: "major cell types (T, B, NK, Monocytes)"
    confidence_threshold: 0.8

  step2:
    method: "reference-based"
    tool: "SingleR"
    reference: "Monaco Immune 2019"
    target: "subtype resolution (CD4/CD8, classical/non-classical monocytes)"
    confidence_threshold: 0.7

  step3:
    method: "llm-assisted"
    tool: "GPT-4"
    target: "ambiguous clusters (clusters 5, 11)"

  validation:
    cross_method_agreement_threshold: 0.9
    manual_review_triggers:
      - "confidence < 0.7"
      - "methods disagree"
      - "cluster size < 50"

# Execution plan
marker_sources:
  - "CellMarker (v2.0)"
  - "PanglaoDB (2024-01)"
  - "User-provided: custom_markers.yaml"

reference_datasets:
  primary: "Monaco_Immune_2019"
  fallback: "HPCA_2018"

estimated_time: "4-6 hours"
confidence_level: "high"

# Risk mitigation
identified_risks:
  - risk: "Cluster 5 may be doublets"
    mitigation: "Flag for manual review, check co-expression of lineage markers"
  - risk: "Rare cell types (<1%) may be missed"
    mitigation: "Use sensitive marker scoring, lower threshold for rare types"

# Agent verdicts
agent_verdicts:
  tissue_biology_expert: "PASS"
  method_specialist: "PASS"
  quality_controller: "CONDITIONAL — address doublet concern for cluster 5"

# Literature added during deliberation
new_literature:
  - "Monaco et al. 2019 — immune cell reference"
  - "Aran et al. 2019 — SingleR method paper"
```

## User Presentation

After Phase 1 completes, present to user:

```markdown
## Phase 1 Complete: Annotation Strategy Approved

**Multi-Agent Deliberation Summary**:
- Rounds: 3
- Outcome: All agents converged on hybrid strategy

**Key Decisions**:
1. **Primary Method**: Hybrid (marker-based + reference-based)
2. **Expected Cell Types**: 8 major types, 15 subtypes
3. **Quality Actions Taken**:
   - Removed 2 low-quality clusters
   - Filtered 3.2% doublets
4. **Estimated Time**: 4-6 hours
5. **Confidence**: HIGH

**Strategy Overview**:
- Step 1: Marker-based for major lineages (T, B, NK, Mono)
- Step 2: SingleR for subtype resolution
- Step 3: LLM-assisted for ambiguous clusters

**Flagged for Manual Review**:
- Cluster 5 (potential doublets)
- Cluster 11 (low marker specificity)

**Next Steps**:
Proceed to Phase 2 (Marker Gene Library Preparation)?

[Approve] [Modify Strategy] [View Detailed Report]
```

## STOP — Wait for User Approval

<IRON-LAW>
After presenting the Phase 1 strategy, **STOP and wait for user approval**.

Do NOT automatically proceed to Phase 2.
Do NOT begin downloading marker databases.
Do NOT start running annotation tools.

The user must explicitly approve the strategy before execution begins.
</IRON-LAW>
