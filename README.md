# MetaBiasSelect
<p align="center">
<img src="https://raw.githubusercontent.com/Sun0shy/MetaBiasSelect/main/man/figures/logo.svg" width="180">
</p>
MetaBiasSelect is an R package for identifying cancer-biased clusters, candidate-biased clusters, and high-confidence shared clusters from multi-cancer single-cell metastatic datasets using a standardized abundance-based workflow.

It is designed for questions such as:

- Which immune clusters are preferentially enriched in one metastatic cancer type?
- Which clusters are shared across multiple metastatic cancers?
- Which clusters should be prioritized for downstream DE, GSEA, trajectory, or communication analysis?

The workflow is **fast** because it mainly operates on **sample × cluster abundance tables**, rather than re-running full gene-level single-cell analysis.

---

# 1. Overview of the workflow

MetaBiasSelect integrates three sources of evidence for each **cluster × cancer** pair:

1. **Cancer-level enrichment**
   - quantified as `log2(O/E)` using sample-level cluster proportions

2. **Cancer-level mean abundance**
   - quantified as the mean cluster proportion within each cancer type

3. **Pairwise differential abundance (DA) support**
   - estimated using pseudobulk `edgeR` across cancer groups **when sufficient biological replication is available**

The package then classifies clusters into:

- **Biased**
- **Candidate_biased**
- **Shared**
- **High-confidence Shared**

---

# 2. Core classification logic

## 2.1 Biased clusters

A cluster is considered **Biased** when:

- DA is available for that cancer, and
- the cluster satisfies at least `min_rules` among:
  - enrichment criterion
  - abundance criterion
  - DA support criterion

Typical default settings:

- `enrich_cutoff = 0.3`
- `prop_cutoff = 0.1`
- `logfc_cutoff = 1`
- `fdr_cutoff = 0.05`
- `min_rules = 2`

---

## 2.2 Candidate-biased clusters

If a cancer group has **insufficient biological replicates**, the pipeline does **not stop**.

Instead, the cluster can still be labeled as **Candidate_biased** when:

- enrichment criterion is met
- abundance criterion is met
- but DA is unavailable because the cancer group has too few samples

This preserves biologically meaningful signals without overclaiming formal statistical support.

---

## 2.3 Shared clusters

A cluster is labeled **Shared** when it is not selected as:

- `Biased`
- `Candidate_biased`

---

## 2.4 High-confidence Shared clusters

A shared cluster is additionally marked as **high-confidence shared** when it satisfies all of the following:

- `classification == "Shared"`
- `shared_score >= shared_score_cutoff`
- `mean_prop >= shared_prop_cutoff`
- present above the abundance cutoff in at least `shared_min_cancers` cancer types

Default settings:

- `shared_score_cutoff = -1.0`
- `shared_prop_cutoff = 0.03`
- `shared_min_cancers = 2`

---

# 3. Heatmap interpretation

The integrated heatmap has:

- one **Shared** column
- one column for each cancer type

## 3.1 Cancer columns

Cancer bias score is based on `log2(O/E)`:

- **red**: enriched in that cancer
- **blue**: depleted in that cancer
- **white**: near baseline

## 3.2 Shared column

Shared stability is based on:

```r
shared_score = -sd(score_across_cancers)
````

Interpretation:

* greener = more stable across cancers
* closer to zero = more shared
* more negative = more variable across cancers

---

# 4. Symbol meanings in the heatmap

The heatmap supports direct symbol annotation.

* `*` = **Biased**
* `†` = **Candidate_biased**
* `•` = **High-confidence Shared**

## Important caution

These symbols do **not** all mean the same thing.

### `*`

This means the cluster is a **formally supported cancer-biased cluster** under the full workflow, including DA support.

### `†`

This means the cluster is a **candidate-biased cluster**:

* enrichment and abundance support it
* but DA was unavailable or insufficient due to limited replication

### `•`

This means the cluster is a **high-confidence shared cluster**:

* not cancer-biased
* stable across cancers
* sufficiently abundant across multiple cancers

So:

* `*` is the strongest evidence for **cancer bias**
* `†` is weaker, but still biologically meaningful
* `•` is **not a significance mark for bias**; it is a mark for **shared stability**

---

# 5. Installation

MetaBiasSelect can be installed either from GitHub or from a local source directory.

## 5.1 Required dependencies

MetaBiasSelect depends on several CRAN packages and one Bioconductor package.

### CRAN packages

- dplyr
- tidyr
- tibble
- ggplot2
- purrr
- rlang
- magrittr
- ggnewscale
- scales
- remotes

### Bioconductor package

- edgeR

---

## 5.2 Install dependencies

Install the required CRAN packages:

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2",
  "purrr", "rlang", "magrittr", "ggnewscale",
  "scales", "remotes"
))
```

Install the required Bioconductor package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")
```

---

## 5.3 Install from GitHub

The recommended installation method is directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Sun0shy/MetaBiasSelect")
```

Then load the package:

```r
library(MetaBiasSelect)
```

---

## 5.4 Install from a local source directory

If you have already downloaded or cloned the repository, you can install it locally.

### From R

```r
remotes::install_local("MetaBiasSelect")
```

or provide the full path:

```r
remotes::install_local("path/to/MetaBiasSelect")
```

### From terminal

```bash
R CMD INSTALL MetaBiasSelect
```

---

## 5.5 Verify installation

After installation, test whether the package can be loaded successfully:

```r
library(MetaBiasSelect)
```

You can also list the main exported functions:

```r
ls("package:MetaBiasSelect")
```

Expected core functions include:

- `run_metabias_pipeline`
- `plot_metabias_heatmap`
- `save_metabias_outputs`
- `run_pairwise_cluster_da`
- `classify_metabias_clusters`

---

## 5.6 Notes

1. `edgeR` must be installed before running the main workflow.
2. If GitHub installation fails because of network or certificate issues, download the repository ZIP file and use `remotes::install_local()` instead.
3. On Linux servers, `R CMD INSTALL MetaBiasSelect` is usually the most stable installation method.
4. After updating the source code, reinstall the package to ensure changes take effect.

---

# 6. Input requirements

MetaBiasSelect accepts either:

## 6.1 A Seurat object

The metadata must contain:

* one column for **sample ID**
* one column for **cancer type**
* one column for **cluster annotation**

Example:

* `sample_id`
* `orig.ident`
* `T_NK_celltype_refined`

## 6.2 A metadata data.frame

Required columns:

* sample column
* cancer column
* cluster column

Each row should represent one cell.

Example:

| sample_id | orig.ident | celltype_annot_final |
| --------- | ---------- | -------------------- |
| GSM1      | CRC        | exhausted CD8 T      |
| GSM1      | CRC        | exhausted CD8 T      |
| GSM2      | PDAC       | memory CD4 T         |

---

# 7. Exported functions and full parameter reference

This package currently exposes five main user-facing functions:

* `run_metabias_pipeline()`
* `plot_metabias_heatmap()`
* `save_metabias_outputs()`
* `run_pairwise_cluster_da()`
* `classify_metabias_clusters()`

---

# 7.1 `run_metabias_pipeline()`

## Purpose

This is the **main workflow function**.

It performs:

1. metadata extraction
2. sample × cluster count table construction
3. sample-level proportion calculation
4. cancer-level enrichment calculation
5. pairwise pseudobulk abundance testing
6. integrated cluster classification
7. score matrix generation
8. optional heatmap generation

## Function signature

```r
run_metabias_pipeline(
  object,
  cluster_col,
  cancer_col,
  sample_col = NULL,
  cancers = NULL,
  enrich_cutoff = 0.3,
  prop_cutoff = 0.1,
  logfc_cutoff = 1,
  fdr_cutoff = 0.05,
  min_rules = 2,
  da_min_support = NULL,
  pseudocount = 1e-6,
  min_samples_per_group = 3,
  cluster_order = NULL,
  return_plot = TRUE,
  shared_score_cutoff = -1.0,
  shared_prop_cutoff = 0.03,
  shared_min_cancers = 2,
  shared_symbol = "•"
)
```

## Parameter explanation

### `object`

Input object.

Can be:

* a **Seurat object**, or
* a **metadata data.frame**

### `cluster_col`

Column name containing cluster annotation.

Examples:

* `T_NK_celltype_refined`
* `celltype_annot_final`

### `cancer_col`

Column name containing cancer type.

Examples:

* `orig.ident`

### `sample_col`

Column name containing sample ID.

If `NULL` or missing, the function will fall back to `cancer_col`, but this is **not recommended** because sample-level replication is crucial for DA analysis.

### `cancers`

Optional vector of cancer types to include and order.

Example:

```r
c("CRC", "PDAC", "STAD")
```

If `NULL`, all cancer labels in the metadata will be used.

### `enrich_cutoff`

Minimum `log2(O/E)` required to support cancer bias.

Higher values are stricter.

### `prop_cutoff`

Minimum cancer-level mean cluster proportion required to support cancer bias.

This prevents very rare clusters from being selected only because of enrichment.

### `logfc_cutoff`

Minimum pairwise DA `logFC` threshold used to count DA support.

Higher values make DA support harder to obtain.

### `fdr_cutoff`

Maximum pairwise DA FDR threshold.

Typical value:

* `0.05`

### `min_rules`

Minimum number of rules passed among:

* enrichment
* proportion
* DA

Typical setting:

* `2`

### `da_min_support`

Minimum number of supporting pairwise DA comparisons.

If `NULL`, all tested other cancers must be outperformed.

Use a numeric value if you want a less strict DA rule.

### `pseudocount`

Small number added to enrichment calculation to avoid division by zero.

Default:

* `1e-6`

Usually no need to change.

### `min_samples_per_group`

Minimum number of samples required in each cancer group for DA testing.

Default:

* `3`

If a cancer has fewer than this number:

* DA is not forced
* the workflow continues
* clusters may still become `Candidate_biased`

This is one of the most important parameters.

### `cluster_order`

Optional manual order for clusters in the final heatmap.

Use this when you want a publication-specific row order.

### `return_plot`

Whether to generate the integrated heatmap.

* `TRUE`: return plot
* `FALSE`: skip plotting

Set `FALSE` if you only want tables or want to debug faster.

### `shared_score_cutoff`

Cutoff for high-confidence shared clusters.

Because `shared_score = -sd(...)`, values closer to zero are more shared.

Typical options:

* `-1.0` = moderate
* `-0.8` = stricter
* `-0.6` = very strict

### `shared_prop_cutoff`

Minimum average abundance required for high-confidence shared clusters.

Typical values:

* `0.03`
* `0.05`

### `shared_min_cancers`

Minimum number of cancer types in which the cluster must exceed the shared abundance threshold.

Typical values:

* `2`
* `3`

### `shared_symbol`

Symbol used on the Shared column for high-confidence shared clusters.

Default:

* `"•"`

---

## Important notes for `run_metabias_pipeline()`

1. **Do not use `sample_col = NULL` unless absolutely necessary.**
   Without real sample IDs, DA analysis is not meaningful.

2. **Small-sample cancers are not dropped.**
   They are downgraded to candidate-only status rather than causing the workflow to fail.

3. **Bias and shared selection are different concepts.**
   A cluster with `•` is not “significant for bias”; it is “stable and shared”.

4. **Thresholds interact.**
   If you tighten `shared_score_cutoff`, fewer shared clusters will get `•`.
   If you tighten `enrich_cutoff` or `prop_cutoff`, fewer biased clusters will get `*`.

---

# 7.2 `plot_metabias_heatmap()`

## Purpose

Draw the integrated heatmap showing:

* cancer bias score
* shared stability
* symbol annotations

## Function signature

```r
plot_metabias_heatmap(
  result,
  cluster_order = NULL,
  class_order = NULL,
  shared_low = "#f7fcf5",
  shared_high = "#31a354",
  bias_low = "#3B4CC0",
  bias_mid = "#F7F7F7",
  bias_high = "#B40426",
  bias_limits = NULL,
  show_symbols = TRUE,
  selected_symbol = "*",
  candidate_symbol = "†",
  shared_symbol = "•",
  symbol_size = 5,
  symbol_color = "black",
  shared_score_cutoff = -1.0,
  shared_prop_cutoff = 0.03,
  shared_min_cancers = 2,
  axis_title_size = 16,
  axis_text_x_size = 13,
  axis_text_y_size = 12,
  legend_title_size = 12,
  legend_text_size = 11,
  axis_line_width = 0.8,
  tile_line_width = 0.8
)
```

## Parameter explanation

### `result`

Output object from `run_metabias_pipeline()`.

### `cluster_order`

Optional manual row order.

### `class_order`

Optional order of classification blocks.

Typical example:

```r
c("Shared", "CRC", "PDAC", "STAD")
```

### `shared_low`, `shared_high`

Low and high colors for the Shared column.

### `bias_low`, `bias_mid`, `bias_high`

Low, middle, and high colors for cancer bias scores.

Default is blue–white–red.

### `bias_limits`

Optional numeric vector of length 2 controlling the displayed range of cancer bias scores.

This is useful when extreme negative values compress the color scale.

Example:

```r
bias_limits = c(-3, 1.5)
```

### `show_symbols`

Whether to overlay symbols.

### `selected_symbol`

Symbol for formally supported biased clusters.

Default:

* `*`

### `candidate_symbol`

Symbol for candidate-biased clusters.

Default:

* `†`

### `shared_symbol`

Symbol for high-confidence shared clusters.

Default:

* `•`

### `symbol_size`

Text size for symbols.

### `symbol_color`

Text color for symbols.

### `shared_score_cutoff`, `shared_prop_cutoff`, `shared_min_cancers`

These determine which Shared clusters get `•`.

These parameters should usually match the values used in `run_metabias_pipeline()`.

### `axis_title_size`

Font size for x/y axis titles.

### `axis_text_x_size`

Font size for x-axis text.

### `axis_text_y_size`

Font size for y-axis text.

### `legend_title_size`

Font size for legend titles.

### `legend_text_size`

Font size for legend text.

### `axis_line_width`

Line width for x/y axes and ticks.

### `tile_line_width`

Border width for heatmap tiles.

---

## Important notes for `plot_metabias_heatmap()`

1. `*`, `†`, and `•` do **not** mean the same type of support.
2. Shared symbols are only placed on the **Shared** column.
3. Candidate-biased clusters should be interpreted more cautiously than fully biased clusters.
4. If your plot looks too compressed, adjust colors, set `bias_limits`, or manually define `cluster_order`.
5. Axis and legend size parameters are intended for publication-ready figure tuning.

---

# 7.3 `save_metabias_outputs()`

## Purpose

Export result tables and plots to disk.

## Function signature

```r
save_metabias_outputs(
  result,
  outdir = "MetaBiasSelect_output",
  prefix = "metabias",
  width = 8,
  height = 6
)
```

## Parameter explanation

### `result`

Output object from `run_metabias_pipeline()`.

### `outdir`

Output directory.

### `prefix`

Prefix added to exported filenames.

### `width`, `height`

Figure width and height for saved heatmaps.

---

## What this function saves

* parameter table
* cancer sample count table
* cluster summary table
* biased cluster table
* shared cluster table
* evidence table
* pairwise DA table
* score matrix
* heatmap PDF
* heatmap PNG
* full RDS object

---

# 7.4 `run_pairwise_cluster_da()`

## Purpose

Run pairwise pseudobulk abundance tests across cancer types using `edgeR`.

## Function signature

```r
run_pairwise_cluster_da(
  count_df,
  sample_meta,
  min_samples_per_group = 3
)
```

## Parameter explanation

### `count_df`

Sample-by-cluster count table.

### `sample_meta`

Sample metadata table with:

* `sample`
* `cancer`

### `min_samples_per_group`

Minimum number of samples required per cancer.

If fewer than two cancer groups meet this requirement, the function returns an empty DA table rather than failing.

---

# 7.5 `classify_metabias_clusters()`

## Purpose

Integrate enrichment, abundance, and DA support into final cluster classification.

## Function signature

```r
classify_metabias_clusters(
  evidence_table,
  enrich_cutoff = 0.3,
  prop_cutoff = 0.1,
  min_rules = 2,
  da_min_support = NULL,
  shared_label = "Shared"
)
```

## Parameter explanation

### `evidence_table`

Cluster × cancer evidence table.

### `enrich_cutoff`

Minimum enrichment threshold.

### `prop_cutoff`

Minimum abundance threshold.

### `min_rules`

Minimum number of rules required for formal bias selection.

### `da_min_support`

Minimum number of supporting pairwise DA contrasts.

### `shared_label`

Label used for non-biased clusters.

---

# 8. Example usage

## 8.1 Seurat object input

```r
library(MetaBiasSelect)

res <- run_metabias_pipeline(
  object = T_cells_merged,
  cluster_col = "T_NK_celltype_refined",
  cancer_col = "orig.ident",
  sample_col = "sample_id",
  cancers = c("CRC", "PDAC", "STAD"),
  enrich_cutoff = 0.3,
  prop_cutoff = 0.10,
  logfc_cutoff = 1,
  fdr_cutoff = 0.05,
  min_rules = 2,
  min_samples_per_group = 3,
  shared_score_cutoff = -1.0,
  shared_prop_cutoff = 0.03,
  shared_min_cancers = 2
)
```

## 8.2 Metadata data.frame input

```r
res <- run_metabias_pipeline(
  object = meta_df,
  cluster_col = "celltype_annot_final",
  cancer_col = "orig.ident",
  sample_col = "sample_id",
  cancers = c("CRC", "PDAC", "STAD")
)
```

---

# 9. How to inspect result tables in R

After running:

```r
res <- run_metabias_pipeline(...)
```

the main output is a list.

## 9.1 See what is inside

```r
names(res)
```

Typical contents:

* `params`
* `meta`
* `sample_meta`
* `cancer_sample_counts`
* `count_df`
* `sample_proportions`
* `enrichment_table`
* `pairwise_da`
* `evidence_table`
* `cluster_summary`
* `biased_clusters_table`
* `shared_clusters_table`
* `score_matrix`
* `plot`

---

## 9.2 View the final summary table

```r
head(res$cluster_summary)
View(res$cluster_summary)
```

Recommended columns to inspect:

```r
res$cluster_summary[, c(
  "cluster",
  "assigned_cancer",
  "classification",
  "selection_status",
  "plot_symbol",
  "log2OE",
  "prop",
  "shared_score",
  "decision_reason"
)]
```

### How to interpret this table

* `classification`: final assigned class
* `selection_status`: Biased / Candidate_biased / Shared
* `plot_symbol`: `*`, `†`, or `•`
* `decision_reason`: why the cluster received that label

---

## 9.3 View the evidence table

```r
head(res$evidence_table)
View(res$evidence_table)
```

Useful columns:

* `Observed`
* `Expected`
* `log2OE`
* `prop`
* `da_support_n`
* `da_tested_n`
* `da_support_prop`
* `cond_enrich`
* `cond_prop`
* `cond_da`
* `selected`
* `candidate_selected`

This is the most important table if you want to verify the evidence behind each cluster label.

---

## 9.4 View pairwise DA results

```r
head(res$pairwise_da)
View(res$pairwise_da)
```

Useful columns:

* `contrast`
* `cancer_high`
* `cancer_low`
* `logFC`
* `FDR`

---

## 9.5 View score matrix

```r
head(res$score_matrix)
View(res$score_matrix)
```

Useful for heatmap input and shared score inspection.

---

## 9.6 View cancer sample counts

```r
res$cancer_sample_counts
```

This table is very important for interpreting `†` symbols.

If a cancer has too few samples, clusters from that cancer may only be eligible for `Candidate_biased`, not full `Biased`.

---

## 9.7 View biased cluster table

```r
head(res$biased_clusters_table)
View(res$biased_clusters_table)
```

This table is useful for directly listing cancer-biased and candidate-biased clusters.

---

## 9.8 View shared cluster table

```r
head(res$shared_clusters_table)
View(res$shared_clusters_table)
```

This table is useful for directly listing high-confidence shared clusters.

---

# 10. How to print and customize the heatmap

## 10.1 Default heatmap

```r
p <- plot_metabias_heatmap(res)
print(p)
```

## 10.2 Customize colors and symbols

```r
p <- plot_metabias_heatmap(
  res,
  shared_low = "#f7fcf5",
  shared_high = "#31a354",
  bias_low = "#3B4CC0",
  bias_mid = "#F7F7F7",
  bias_high = "#B40426",
  selected_symbol = "*",
  candidate_symbol = "†",
  shared_symbol = "•",
  symbol_size = 5,
  symbol_color = "black"
)

print(p)
```

## 10.3 Use publication-oriented bias scaling

```r
p <- plot_metabias_heatmap(
  res,
  bias_limits = c(-3, 1.5)
)

print(p)
```

This is useful when extreme negative values dominate the color scale.

---

# 11. How to save results

Use:

```r
save_metabias_outputs(
  result = res,
  outdir = "MetaBiasSelect_output",
  prefix = "T_NK"
)
```

This will create files such as:

* `T_NK_parameters.csv`
* `T_NK_cancer_sample_counts.csv`
* `T_NK_cluster_summary.csv`
* `T_NK_biased_clusters_table.csv`
* `T_NK_shared_clusters_table.csv`
* `T_NK_cluster_cancer_evidence.csv`
* `T_NK_pairwise_DA.csv`
* `T_NK_score_matrix.csv`
* `T_NK_heatmap.pdf`
* `T_NK_heatmap.png`
* `T_NK_full_result.rds`

---

# 12. Recommended interpretation strategy

## 12.1 Highest-confidence cancer-biased clusters

Use clusters marked with `*`.

These are the strongest candidates for downstream biological interpretation.

## 12.2 Candidate-biased clusters

Use clusters marked with `†`.

These are useful when:

* the biological pattern is clear
* but replication is insufficient for full DA support

Interpret them cautiously.

## 12.3 Highest-confidence shared clusters

Use clusters marked with `•`.

These are the best candidates for the **shared metastatic immune background**.

---

# 13. Common pitfalls and cautions

## 13.1 Using cancer ID as sample ID

Do not do this unless you have no alternative.
Without true sample-level replication, DA loses meaning.

## 13.2 Over-interpreting `†`

`†` does not mean full statistical support.
It means “biologically plausible candidate under limited replication”.

## 13.3 Treating `•` like a p-value mark

Do not do this.
`•` is a **shared stability mark**, not a bias significance mark.

## 13.4 Choosing thresholds too aggressively

If thresholds are too strict:

* almost no clusters will get symbols

If thresholds are too loose:

* many weak clusters will be selected

A practical strategy is:

* start with defaults
* inspect `cluster_summary`
* then tune according to dataset size and biology

## 13.5 Forgetting to check sample counts

Always inspect:

```r
res$cancer_sample_counts
```

This tells you whether `*` is even possible for each cancer.

---

# 14. Minimal reproducible example

```r
library(MetaBiasSelect)

res <- run_metabias_pipeline(
  object = T_cells_merged,
  cluster_col = "T_NK_celltype_refined",
  cancer_col = "orig.ident",
  sample_col = "sample_id",
  cancers = c("CRC", "PDAC", "STAD")
)

print(res$plot)

head(res$cluster_summary)
head(res$evidence_table)

save_metabias_outputs(
  result = res,
  outdir = "MetaBiasSelect_output",
  prefix = "example"
)
```

---

# 15. Suggested reporting language

If you describe this workflow in a manuscript, a suitable wording is:

**We used a standardized multi-cancer metastatic bias-cluster selection workflow integrating cancer-level enrichment, cancer-level mean abundance, and pairwise differential abundance support to identify cancer-biased, candidate-biased, and shared cell states.**

---

# 16. License

This software is currently provided for academic research use.

For redistribution, modification, or commercial use, please contact the author.

---

# 17. Author

Hongyu Sun
