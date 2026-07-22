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

## Examples

``` r
# use a sequencing result example
seq_results <- example("Sequencing results")

# plotting the depth ratio over all samples
plot_DR(seq_results)


# plotting the depth ratio for sample S_2_2 only
plot_DR(seq_results, samples = "S_2_2")


# plotting the depth ratio for S_2_2 with labels
plot_DR(seq_results, samples = "S_2_2",
        labels = seq_results$mutations["nature"])


# let us define a function to filter germinal and pre-neoplastic
# from the input data
library(dplyr)
filter_data <- function(data) {
  data %>% dplyr::filter(!nature %in% list("germinal",
                                           "pre-neoplastic"))
}

# plotting the depth ratio without germinal and pre-neoplastic
plot_DR(seq_results, mutation_filter = filter_data)


# the same plots can be drawn by using the mutations data frame
# in place of the `simulate_seq()` output

# filter germinal mutations
f_seq <- seq_results$mutations %>%
   dplyr::filter(nature != "germinal")

# plotting the depth ratio
plot_DR(f_seq)

```
