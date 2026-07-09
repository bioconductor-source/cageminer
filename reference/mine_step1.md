# Step 1: Get all putative candidate genes for a given sliding window

For a user-defined sliding window relative to each SNP, this function
will subset all genes whose genomic positions overlap with the sliding
window.

## Usage

``` r
mine_step1(gene_ranges, marker_ranges, window = 2, expand_intervals = TRUE)
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

## Value

A GRanges or GRangesList object with the genomic positions of all
putative candidate genes.

## See also

[`findOverlaps-methods`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)

## Examples

``` r
data(snp_pos)
data(gene_ranges)
genes <- mine_step1(gene_ranges, snp_pos, window = 2)
```
