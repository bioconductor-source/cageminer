# Plot a barplot of SNP distribution across chromosomes

Plot a barplot of SNP distribution across chromosomes

## Usage

``` r
plot_snp_distribution(marker_ranges)
```

## Arguments

- marker_ranges:

  Genomic positions of SNPs. For a single trait, a GRanges object. For
  multiple traits, a GRangesList or CompressedGRangesList object, with
  each element of the list representing SNP positions for a particular
  trait. List elements must have names for proper labelling.

## Value

A ggplot object.

## Examples

``` r
data(snp_pos)
p <- plot_snp_distribution(snp_pos)
```
