# Step 3: Select candidates based on gene significance

Step 3: Select candidates based on gene significance

## Usage

``` r
mine_step3(
  exp,
  metadata,
  metadata_cols = 1,
  candidates,
  sample_group,
  min_cor = 0.2,
  alpha = 0.05,
  ...
)
```

## Arguments

- exp:

  Expression data frame with genes in row names and samples in column
  names or a SummarizedExperiment object.

- metadata:

  Sample metadata with samples in row names and sample information in
  the first column. Ignored if `exp` is a SummarizedExperiment object,
  as the colData will be extracted from the object.

- metadata_cols:

  A vector (either numeric or character) indicating which columns should
  be extracted from column metadata if **exp** is a
  `SummarizedExperiment` object. The vector can contain column indices
  (numeric) or column names (character). By default, all columns are
  used.

- candidates:

  Character vector of candidate genes to be inspected.

- sample_group:

  Level of sample metadata to be used for filtering in gene-trait
  correlation.

- min_cor:

  Minimum correlation value for
  [`BioNERO::gene_significance()`](https://rdrr.io/pkg/BioNERO/man/gene_significance.html).
  Default: 0.2

- alpha:

  Numeric indicating significance level. Default: 0.05

- ...:

  Additional arguments to
  [`BioNERO::gene_significance`](https://rdrr.io/pkg/BioNERO/man/gene_significance.html).

## Value

A data frame with mined candidate genes and their correlation to the
condition of interest.

## Examples

``` r
data(pepper_se)
data(snp_pos)
data(gene_ranges)
data(guides)
data(gcn)
data(mine2)
set.seed(1)
mine3 <- mine_step3(
    exp = pepper_se,
    candidates = mine2$candidates,
    sample_group = "PRR_stress"
)
```
