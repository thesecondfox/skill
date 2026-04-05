---
name: bio-pathway-kegg-pathways
description: KEGG pathway and module enrichment analysis using clusterProfiler enrichKEGG and enrichMKEGG. Use when identifying metabolic and signaling pathways over-represented in a gene list. Supports 4000+ organisms via KEGG online database. Supports local-first mode via local-reference-registry.
tool_type: r
primary_tool: clusterProfiler
depends_on: local-reference-registry
---

## Version Compatibility

Reference examples tested with: R stats (base), clusterProfiler 4.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# KEGG Pathway Enrichment

## Local-First Mode (PRIORITY)

**Before running any KEGG enrichment, check for local reference data:**

```r
config_path <- path.expand("~/.bioref_config.yaml")

if (file.exists(config_path)) {
    # LOCAL MODE: Use pre-downloaded KEGG GMT files (no network needed)
    library(clusterProfiler)
    library(yaml)

    config <- read_yaml(config_path)
    ref_root <- config$ref_root

    # Load local KEGG gene sets
    gmt_kegg <- file.path(ref_root, config$gmt_files$human$kegg)  # or $mouse$kegg
    kegg_terms <- read.gmt(gmt_kegg)

    # ORA with local KEGG (uses enricher instead of enrichKEGG)
    kk <- enricher(
        gene = gene_list,            # Gene symbols (matching GMT format)
        TERM2GENE = kegg_terms,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
    )

    # GSEA with local KEGG
    gsea_kegg <- GSEA(
        geneList = ranked_gene_list,
        TERM2GENE = kegg_terms,
        pvalueCutoff = 0.05
    )

} else {
    # ONLINE MODE: Fallback to enrichKEGG (requires network)
    message("No local reference config (~/.bioref_config.yaml) found.")
    message("Using online mode. If network is unavailable, run local-reference-registry first.")
}
```

> **Note:** Local mode uses `enricher()` / `GSEA()` + `TERM2GENE` instead of `enrichKEGG()` / `gseKEGG()`. This eliminates network dependency. Gene IDs in GMT files are typically symbols.

## Online Mode (Original Pattern)

## Core Pattern

**Goal:** Identify KEGG metabolic and signaling pathways over-represented in a gene list.

**Approach:** Test for enrichment using the hypergeometric test via clusterProfiler enrichKEGG against the KEGG online database.

**"Find enriched KEGG pathways in my gene list"** → Test whether KEGG pathway gene sets are over-represented among significant genes.

```r
library(clusterProfiler)

kk <- enrichKEGG(
    gene = gene_list,           # Character vector of gene IDs
    organism = 'hsa',           # KEGG organism code
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)
```

## Prepare Gene List

**Goal:** Extract significant Entrez gene IDs from DE results in the format required by enrichKEGG.

**Approach:** Filter by significance thresholds and convert gene symbols to Entrez IDs (KEGG requires NCBI Entrez).

```r
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

# KEGG requires NCBI Entrez gene IDs (kegg, ncbi-geneid)
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

## KEGG ID Conversion

**Goal:** Convert between KEGG-specific identifiers and other gene ID formats.

**Approach:** Use bitr_kegg to map between kegg, ncbi-geneid, ncbi-proteinid, and uniprot ID types.

```r
# Convert between KEGG and other IDs
kegg_ids <- bitr_kegg(gene_list, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')

# Available types: kegg, ncbi-geneid, ncbi-proteinid, uniprot
```

## Run KEGG Pathway Enrichment

**Goal:** Perform KEGG pathway over-representation analysis with customizable parameters.

**Approach:** Run enrichKEGG with specified organism, ID type, and statistical thresholds.

```r
kk <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',
    keyType = 'ncbi-geneid',    # or 'kegg'
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = 500
)

# View results
head(kk)
results <- as.data.frame(kk)
```

## Make Results Readable

```r
# enrichKEGG does NOT have readable parameter - use setReadable
library(org.Hs.eg.db)
kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```

## KEGG Module Enrichment

**Goal:** Test for enrichment of KEGG modules (smaller functional units than pathways).

**Approach:** Use enrichMKEGG which tests against KEGG module definitions rather than full pathways.

```r
# KEGG modules are smaller functional units than pathways
mkk <- enrichMKEGG(
    gene = gene_list,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Common Organism Codes

| Organism | Code | Common Name |
|----------|------|-------------|
| hsa | Human | Homo sapiens |
| mmu | Mouse | Mus musculus |
| rno | Rat | Rattus norvegicus |
| dre | Zebrafish | Danio rerio |
| dme | Fruit fly | Drosophila melanogaster |
| cel | Worm | C. elegans |
| sce | Yeast | S. cerevisiae |
| ath | Arabidopsis | A. thaliana |
| eco | E. coli K-12 | |

```r
# Find organism codes
search_kegg_organism('mouse')
search_kegg_organism('zebrafish')
```

## With Background Universe

**Goal:** Restrict KEGG enrichment to genes actually measured in the experiment.

**Approach:** Convert all tested genes to Entrez IDs and pass as the universe parameter.

```r
all_genes <- de_results$gene_id
universe_ids <- bitr(all_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(
    gene = gene_list,
    universe = universe_ids$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Extract and Export Results

**Goal:** Save KEGG enrichment results to CSV and extract genes belonging to specific pathways.

**Approach:** Convert enrichment object to data frame, export, and access pathway gene sets via the geneSets slot.

```r
# Convert to data frame
results_df <- as.data.frame(kk)

# Key columns: ID (pathway), Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID, Count

# Export
write.csv(results_df, 'kegg_enrichment_results.csv', row.names = FALSE)

# Get genes in a specific pathway
pathway_genes <- kk@geneSets[['hsa04110']]  # Cell cycle
```

## Browse KEGG Pathways

**Goal:** Visualize enriched genes overlaid on KEGG pathway diagrams.

**Approach:** Use browseKEGG for interactive browser view or pathview to generate annotated pathway images.

```r
# View pathway in browser (opens KEGG website)
browseKEGG(kk, 'hsa04110')

# Download pathway image
library(pathview)
pathview(gene.data = gene_list, pathway.id = 'hsa04110', species = 'hsa')
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| organism | hsa | KEGG organism code |
| keyType | kegg | Input ID type |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value threshold |
| pAdjustMethod | BH | Adjustment method |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per pathway |
| maxGSSize | 500 | Max genes per pathway |
| use_internal_data | FALSE | Use local KEGG data |

## Compare Multiple Gene Lists

**Goal:** Compare KEGG pathway enrichment across multiple gene lists (e.g., upregulated vs downregulated).

**Approach:** Use compareCluster with enrichKEGG to run enrichment per group and visualize with dotplot.

```r
# Compare KEGG enrichment across groups
gene_lists <- list(
    up = up_genes,
    down = down_genes
)

ck <- compareCluster(
    geneClusters = gene_lists,
    fun = 'enrichKEGG',
    organism = 'hsa'
)

dotplot(ck)
```

## Notes

- **No readable parameter** - use `setReadable()` with OrgDb
- **Requires internet** - queries KEGG database online
- **use_internal_data** - set TRUE to use cached KEGG data (may be outdated)
- **Pathway IDs** - format is organism code + 5 digits (e.g., hsa04110)

## Related Skills

- go-enrichment - Gene Ontology enrichment analysis
- gsea - GSEA using KEGG pathways (gseKEGG)
- enrichment-visualization - Visualize KEGG results
