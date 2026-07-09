# Circos plot of SNP distribution across chromosomes

Circos plot of SNP distribution across chromosomes

## Usage

``` r
plot_snp_circos(genome_ranges, gene_ranges, marker_ranges)
```

## Arguments

- genome_ranges:

  A GRanges object with chromosome lengths.

- gene_ranges:

  A GRanges object with genomic coordinates of all genes in the genome.

- marker_ranges:

  Genomic positions of SNPs. For a single trait, a GRanges object. For
  multiple traits, a GRangesList or CompressedGRangesList object, with
  each element of the list representing SNP positions for a particular
  trait.

## Value

A ggplot object with a circos plot of molecular marker distribution
across chromosomes.

## Examples

``` r
data(snp_pos)
data(gene_ranges)
data(chr_length)
p <- plot_snp_circos(chr_length, gene_ranges, snp_pos)
```
