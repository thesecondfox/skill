---
name: local-reference-registry
description: Local reference data registry for offline bioinformatics analysis. Manages GMT gene sets, organism annotation databases, and ligand-receptor databases. Eliminates network dependency for clusterProfiler, GSEA, CellChat and other tools by providing a unified local data path configuration.
tool_type: universal
primary_tool: configuration
---

# Local Reference Registry

## Persona and Role

You are a **Reference Data Manager** responsible for configuring and maintaining local bioinformatics reference data paths. Your primary goal is to eliminate network dependencies by directing all downstream analysis skills to use locally stored gene sets, annotation databases, and other reference files.

## When This Skill Activates

This skill MUST be checked **before** any of the following operations:
- GO enrichment analysis (enrichGO, gseGO)
- KEGG pathway analysis (enrichKEGG, gseKEGG)
- GSEA with MSigDB gene sets
- Reactome pathway analysis
- CellChat / NicheNet cell communication analysis
- Any analysis requiring GMT files or organism annotation databases

## Configuration File

Location: `~/.bioref_config.yaml`

This file is the single source of truth for all local reference data paths. Any skill that needs reference data reads this file first.

### File Format

```yaml
# ~/.bioref_config.yaml
# Local Reference Data Registry
# Generated or manually edited

ref_root: /path/to/your/reference/data

# GMT gene set files (MSigDB format)
gmt_files:
  human:
    hallmark: h.all.v2024.1.Hs.symbols.gmt
    kegg: c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt
    reactome: c2.cp.reactome.v2024.1.Hs.symbols.gmt
    go_bp: c5.go.bp.v2024.1.Hs.symbols.gmt
    go_mf: c5.go.mf.v2024.1.Hs.symbols.gmt
    go_cc: c5.go.cc.v2024.1.Hs.symbols.gmt
    wikipathways: c2.cp.wikipathways.v2024.1.Hs.symbols.gmt
    tft: c3.tft.v2024.1.Hs.symbols.gmt
    immunesigdb: c7.immunesigdb.v2024.1.Hs.symbols.gmt
    oncogenic: c6.all.v2024.1.Hs.symbols.gmt
    positional: c1.all.v2024.1.Hs.symbols.gmt
  mouse:
    hallmark: mh.all.v2024.1.Mm.symbols.gmt
    kegg: m2.cp.kegg_medicus.v2024.1.Mm.symbols.gmt
    reactome: m2.cp.reactome.v2024.1.Mm.symbols.gmt
    go_bp: m5.go.bp.v2024.1.Mm.symbols.gmt
    go_mf: m5.go.mf.v2024.1.Mm.symbols.gmt
    go_cc: m5.go.cc.v2024.1.Mm.symbols.gmt

# Organism annotation databases (OrgDb)
orgdb:
  human: org.Hs.eg.db
  mouse: org.Mm.eg.db

# CellChat ligand-receptor databases
cellchat_db:
  human: CellChatDB.human.rds
  mouse: CellChatDB.mouse.rds

# Custom gene sets (user-defined)
custom_gmt: []
```

## Workflow

### Step 1: Check Configuration

Before any enrichment or pathway analysis, check if the config file exists:

```python
import os
import yaml

config_path = os.path.expanduser("~/.bioref_config.yaml")

if os.path.exists(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)
    ref_root = config.get("ref_root", "")
    print(f"Local reference data found: {ref_root}")
else:
    print("No local reference config found.")
    print("Please provide your local reference data path,")
    print("or I will use online mode (requires network).")
```

```r
config_path <- path.expand("~/.bioref_config.yaml")

if (file.exists(config_path)) {
  config <- yaml::read_yaml(config_path)
  ref_root <- config$ref_root
  message("Local reference data found: ", ref_root)
} else {
  message("No local reference config found.")
  message("Please provide your local reference data path,")
  message("or I will use online mode (requires network).")
}
```

### Step 2: Initialize Configuration (First Time)

When a user provides a reference data path for the first time:

```python
import os
import glob
import yaml

def init_reference_registry(ref_root):
    """Scan reference directory and generate config."""
    ref_root = os.path.expanduser(ref_root)

    if not os.path.isdir(ref_root):
        raise FileNotFoundError(f"Path does not exist: {ref_root}")

    # Scan for GMT files
    gmt_files = glob.glob(os.path.join(ref_root, "**/*.gmt"), recursive=True)

    config = {
        "ref_root": ref_root,
        "gmt_files": {"human": {}, "mouse": {}},
        "orgdb": {"human": "org.Hs.eg.db", "mouse": "org.Mm.eg.db"},
        "cellchat_db": {"human": "", "mouse": ""},
        "custom_gmt": []
    }

    # Auto-classify GMT files
    for gmt in gmt_files:
        basename = os.path.basename(gmt)
        relpath = os.path.relpath(gmt, ref_root)

        # Determine species
        if ".Hs." in basename or "human" in basename.lower():
            species = "human"
        elif ".Mm." in basename or "mouse" in basename.lower():
            species = "mouse"
        else:
            config["custom_gmt"].append(relpath)
            continue

        # Determine category
        basename_lower = basename.lower()
        if "hallmark" in basename_lower or basename.startswith("h.all"):
            config["gmt_files"][species]["hallmark"] = relpath
        elif "kegg" in basename_lower:
            config["gmt_files"][species]["kegg"] = relpath
        elif "reactome" in basename_lower:
            config["gmt_files"][species]["reactome"] = relpath
        elif "go.bp" in basename_lower or "gobp" in basename_lower:
            config["gmt_files"][species]["go_bp"] = relpath
        elif "go.mf" in basename_lower or "gomf" in basename_lower:
            config["gmt_files"][species]["go_mf"] = relpath
        elif "go.cc" in basename_lower or "gocc" in basename_lower:
            config["gmt_files"][species]["go_cc"] = relpath
        elif "wikipathways" in basename_lower:
            config["gmt_files"][species]["wikipathways"] = relpath
        elif "tft" in basename_lower:
            config["gmt_files"][species]["tft"] = relpath
        elif "immunesigdb" in basename_lower:
            config["gmt_files"][species]["immunesigdb"] = relpath
        elif "oncogenic" in basename_lower or "c6" in basename_lower:
            config["gmt_files"][species]["oncogenic"] = relpath
        elif "positional" in basename_lower or "c1" in basename_lower:
            config["gmt_files"][species]["positional"] = relpath
        else:
            config["custom_gmt"].append(relpath)

    # Scan for CellChat DB
    for rds in glob.glob(os.path.join(ref_root, "**/*.rds"), recursive=True):
        rds_name = os.path.basename(rds).lower()
        relpath = os.path.relpath(rds, ref_root)
        if "cellchat" in rds_name and "human" in rds_name:
            config["cellchat_db"]["human"] = relpath
        elif "cellchat" in rds_name and "mouse" in rds_name:
            config["cellchat_db"]["mouse"] = relpath

    # Write config
    config_path = os.path.expanduser("~/.bioref_config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)

    print(f"Config saved to {config_path}")
    print(f"Found {len(gmt_files)} GMT files")
    print(f"Human gene sets: {len(config['gmt_files']['human'])}")
    print(f"Mouse gene sets: {len(config['gmt_files']['mouse'])}")
    print(f"Custom/unclassified: {len(config['custom_gmt'])}")

    return config
```

### Step 3: Load GMT in Downstream Skills

#### R (clusterProfiler)

```r
# Load local GMT instead of using enrichGO/enrichKEGG online
library(clusterProfiler)
library(yaml)

config <- read_yaml(path.expand("~/.bioref_config.yaml"))
ref_root <- config$ref_root

# --- GO Enrichment (local mode) ---
gmt_go_bp <- file.path(ref_root, config$gmt_files$human$go_bp)
go_terms <- read.gmt(gmt_go_bp)

ego <- enricher(
    gene = gene_list,
    TERM2GENE = go_terms,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

# --- KEGG Enrichment (local mode) ---
gmt_kegg <- file.path(ref_root, config$gmt_files$human$kegg)
kegg_terms <- read.gmt(gmt_kegg)

ekegg <- enricher(
    gene = gene_list,
    TERM2GENE = kegg_terms,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)

# --- GSEA (local mode) ---
gmt_hallmark <- file.path(ref_root, config$gmt_files$human$hallmark)
hallmark_terms <- read.gmt(gmt_hallmark)

gsea_result <- GSEA(
    geneList = ranked_gene_list,
    TERM2GENE = hallmark_terms,
    pvalueCutoff = 0.05
)
```

#### Python (gseapy)

```python
import gseapy as gp
import yaml
import os

with open(os.path.expanduser("~/.bioref_config.yaml")) as f:
    config = yaml.safe_load(f)

ref_root = config["ref_root"]

# GSEA with local GMT
gmt_path = os.path.join(ref_root, config["gmt_files"]["human"]["hallmark"])

gsea_results = gp.gsea(
    data=expression_df,
    gene_sets=gmt_path,
    permutation_num=1000,
    outdir="gsea_output"
)

# ORA with local GMT
ora_results = gp.enrich(
    gene_list=gene_list,
    gene_sets=gmt_path,
    outdir="ora_output"
)
```

### Step 4: CellChat Local Database

```r
library(yaml)

config <- read_yaml(path.expand("~/.bioref_config.yaml"))
ref_root <- config$ref_root

# Load local CellChat database
cellchat_db_path <- file.path(ref_root, config$cellchat_db$human)
if (file.exists(cellchat_db_path)) {
    CellChatDB <- readRDS(cellchat_db_path)
    message("Using local CellChat database")
} else {
    CellChatDB <- CellChatDB.human
    message("Local DB not found, using built-in CellChat database")
}

cellchat@DB <- CellChatDB
```

## Behavior Rules for All Downstream Skills

1. **Always check `~/.bioref_config.yaml` first** before any enrichment, pathway, or gene set analysis
2. **If config exists**: use local GMT files via `enricher()` / `GSEA()` with `TERM2GENE` parameter instead of `enrichGO()` / `enrichKEGG()` which require network
3. **If config does not exist**: ask the user once — "Do you have local reference data (GMT files)? If yes, provide the path. If no, I will use online mode."
4. **Never hardcode paths**: always read from config
5. **Species awareness**: check whether the analysis is human or mouse, load the corresponding GMT files
6. **Fallback gracefully**: if a specific GMT category is missing from config, warn the user and offer online fallback

## Supported GMT Categories

| Category Key | MSigDB Collection | Description |
|-------------|-------------------|-------------|
| hallmark | H | Hallmark gene sets (50 sets) |
| kegg | C2:KEGG | KEGG pathway gene sets |
| reactome | C2:Reactome | Reactome pathway gene sets |
| go_bp | C5:GO:BP | Gene Ontology Biological Process |
| go_mf | C5:GO:MF | Gene Ontology Molecular Function |
| go_cc | C5:GO:CC | Gene Ontology Cellular Component |
| wikipathways | C2:WikiPathways | WikiPathways gene sets |
| tft | C3:TFT | Transcription factor targets |
| immunesigdb | C7:ImmuneSigDB | Immunologic signature gene sets |
| oncogenic | C6 | Oncogenic signature gene sets |
| positional | C1 | Positional gene sets |

## Quick Reference

```bash
# Initialize config (run once)
python -c "
from local_reference_registry import init_reference_registry
init_reference_registry('/your/path/to/Analysis_ref')
"

# Verify config
cat ~/.bioref_config.yaml

# Update path
python -c "
import yaml
config = yaml.safe_load(open(os.path.expanduser('~/.bioref_config.yaml')))
config['ref_root'] = '/new/path'
yaml.dump(config, open(os.path.expanduser('~/.bioref_config.yaml'), 'w'))
"
```
