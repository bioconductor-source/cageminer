# Score candidate genes and select the top n genes

Score candidate genes and select the top n genes

## Usage

``` r
score_genes(
  mined_candidates,
  hubs = NULL,
  tfs = NULL,
  pick_top = 10,
  weight_tf = 2,
  weight_hub = 2,
  weight_both = 3
)
```

## Arguments

- mined_candidates:

  Data frame resulting from [`mine_candidates()`](mine_candidates.md) or
  `mine_step()`.

- hubs:

  Character vector of hub genes.

- tfs:

  Character vector of transcription factors.

- pick_top:

  Number of top genes to select. Default: 10.

- weight_tf:

  Numeric scalar with the weight to which correlation coefficients will
  be multiplied if the gene is a TF. Default: 2.

- weight_hub:

  Numeric scalar with the weight to which correlation coefficients will
  be multiplied if the gene is a hub. Default: 2.

- weight_both:

  Numeric scalar with the weight to which correlation coefficients will
  be multiplied if the gene is both a TF and a hub. Default: 3.

## Value

Data frame with top n candidates and their scores.

## Examples

``` r
data(tfs)
data(hubs)
data(mined_candidates)
set.seed(1)
scored <- score_genes(mined_candidates, hubs$Gene, tfs$Gene_ID)
#> Number of genes < 'pick_top'. Picking all genes.
```
