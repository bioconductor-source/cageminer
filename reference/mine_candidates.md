# Mine high-confidence candidate genes in a single step

Mine high-confidence candidate genes in a single step

## Usage

``` r
mine_candidates(
  gene_ranges = NULL,
  marker_ranges = NULL,
  window = 2,
  expand_intervals = TRUE,
  gene_col = "ID",
  exp = NULL,
  gcn = NULL,
  guides = NULL,
  metadata,
  metadata_cols = 1,
  sample_group,
  min_cor = 0.2,
  alpha = 0.05,
  ...
)
```

## Arguments

- gene_ranges:

  A GRanges object with genomic coordinates of all genes in the genome.

- marker_ranges:

  Genomic positions of SNPs. For a single trait, a GRanges object. For
  multiple traits, a GRangesList or CompressedGRangesList object, with
  each element of the list representing SNP positions for a particular
  trait.

- window:

  Sliding window (in Mb) upstream and downstream relative to each SNP.
  Default: 2.

- expand_intervals:

  Logical indicating whether or not to expand markers that are
  represented by intervals. This is particularly useful if users want to
  use a custom interval defined by linkage disequilibrium, for example.
  Default: TRUE.

- gene_col:

  Column of the GRanges object containing gene ID. Default: "ID", the
  default for gff/gff3 files imported with rtracklayer::import.

- exp:

  Expression data frame with genes in row names and samples in column
  names or a SummarizedExperiment object.

- gcn:

  Gene coexpression network returned by
  [`BioNERO::exp2gcn()`](https://rdrr.io/pkg/BioNERO/man/exp2gcn.html).

- guides:

  Guide genes as a character vector or as a data frame with genes in the
  first column and gene annotation class in the second column.

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
# \donttest{
data(pepper_se)
data(snp_pos)
data(gene_ranges)
data(guides)
data(gcn)
set.seed(1)
candidates <- mine_candidates(gene_ranges, snp_pos, exp = pepper_se,
                              gcn = gcn, guides = guides$Gene,
                              sample_group = "PRR_stress")
#> Enrichment analysis for module cyan...
# }
```
