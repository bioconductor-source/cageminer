# Step 2: Get candidates in modules enriched in guide genes

Step 2: Get candidates in modules enriched in guide genes

## Usage

``` r
mine_step2(exp, gcn, guides, candidates, ...)
```

## Arguments

- exp:

  Expression data frame with genes in row names and samples in column
  names or a SummarizedExperiment object.

- gcn:

  Gene coexpression network returned by
  [`BioNERO::exp2gcn()`](https://rdrr.io/pkg/BioNERO/man/exp2gcn.html).

- guides:

  Guide genes as a character vector or as a data frame with genes in the
  first column and gene annotation class in the second column.

- candidates:

  Character vector of all candidates genes to be inspected.

- ...:

  Additional arguments to
  [`BioNERO::module_enrichment`](https://rdrr.io/pkg/BioNERO/man/module_enrichment.html)

## Value

A list of 2 elements:

- candidates:

  Character vector of candidates after step 2

- enrichment:

  Data frame of results for enrichment analysis

## Examples

``` r
data(pepper_se)
data(guides)
data(gcn)
set.seed(1)
mine2 <- mine_step2(
    exp = pepper_se,
    gcn = gcn,
    guides = guides$Gene,
    candidates = rownames(pepper_se)
)
#> Enrichment analysis for module cyan...
```
