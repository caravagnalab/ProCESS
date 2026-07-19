# Plot the genome-wide Depth Ratio (DR)

This function plots the genome-wide Depth Ratio (DR) of a specific
sample.

## Usage

``` r
plot_DR(
  seq_result,
  chromosomes = NULL,
  samples = NULL,
  labels = NULL,
  mutation_filter = filter_germinal,
  cuts = c(0, 1),
  N = 5000
)
```

## Arguments

- seq_result:

  Either the output of
  [`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md)
  or a data frame containing sequencing results.

- chromosomes:

  A character vector specifying the chromosomes to include in the plot
  (default: all the chromosomes in `seq_res`).

- samples:

  A character vector specifying the sample names to include in the plot.
  When set to `NULL`, the function includes all samples except the
  "normal sample" (default: `NULL`).

- labels:

  A data frame column labelling each mutation. Each label is associated
  to a different colour in the plot (default: `NULL`).

- mutation_filter:

  A function filtering mutations from the input data (default: a
  function filtering out "germinal" mutations, e.g.,
  `function(x) x %>% dplyr::filter(nature != "germinal")`).

- cuts:

  A numeric vector specifying the range of VAF values to include in the
  plot (default: `c(0, 1)`).

- N:

  The number of mutations to sample for plotting (default: 5000).

## Value

A ggplot2 object showing the DR distribution across the genome.

## See also

[`plot_VAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF.md),
[`plot_BAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_BAF.md)
