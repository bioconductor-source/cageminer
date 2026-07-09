# Simulation of the output list from BioNERO::exp2gcn() with pepper data

This object is a list as returned by BioNERO::exp2gcn(), but only the
element genes_and_modules is included. For running time issues, only
genes in the cyan module were kept in the element genes_and_modules. All
other list elements have been assigned NULL. The network was inferred
using the code from the vignette.

## Usage

``` r
data(gcn)
```

## Format

A list with the elements returned by BioNERO::exp2gcn().

## Examples

``` r
data(gcn)
```
