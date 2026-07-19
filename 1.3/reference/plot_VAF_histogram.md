# Plot a Variant Allele Frequency (VAF) histogram

This function generates a histogram showing the distribution of Variant
Allele Frequency (VAF) across samples and chromosomes.

## Usage

``` r
plot_VAF_histogram(
  seq_result,
  chromosomes = NULL,
  samples = NULL,
  labels = NULL,
  binwidth = NULL,
  mutation_filter = filter_germinal,
  driver_mutations = NULL,
  driver_mutation_labels = TRUE,
  cuts = c(0, 1)
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

- binwidth:

  The width of the plot bins. When set to `NULL`, the function computes
  the most convenient bin width according to the maximum coverage
  reported in the data frame (default: `NULL`).

- mutation_filter:

  A function filtering mutations from the input data (default: a
  function filtering out "germinal" mutations, e.g.,
  `function(x) x %>% dplyr::filter(nature != "germinal")`).

- driver_mutations:

  The data frame of the driver mutations as returned by
  `PhylogeneticForest$get_driver_mutations()`. This parameter can be
  avoided when `seq_result` is the result of the function
  [`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md)
  (default: `NULL`).

- driver_mutation_labels:

  A Boolean value to enable/disable driver mutation labels in the
  returned plot (default: TRUE).

- cuts:

  A numeric vector specifying the range of VAF values to include in the
  plot (default: `c(0, 1)`).

## Value

A ggplot2 object showing the VAF histogram.

## See also

[`plot_VAF_marginals()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF_marginals.md),
[`plot_VAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF.md)

## Examples

``` r
# use a sequencing result example
seq_results <- example("Sequencing results")

# plotting the VAF histogram without germinal mutations
plot_VAF_histogram(seq_results)


# let us define a function to filter germinal and pre-neoplastic
# from the input data
library(dplyr)
filter_data <- function(data) {
  data %>% dplyr::filter(!nature %in% list("germinal",
                                           "pre-neoplastic"))
}

# plotting the VAF histogram without germinal and pre-neoplastic
plot_VAF_histogram(seq_results, mutation_filter = filter_data)


# plotting the VAF histogram filtering out VAFs below 0.02
plot_VAF_histogram(seq_results, cuts = c(0.02, 1))


# plotting the VAF histogram with labels
plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
                   labels = seq_results$mutations["cause"])


# avoid the driver mutation labels
plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
                   labels = seq_results$mutations["cause"],
                   driver_mutation_labels = FALSE)


# the same plots can be drawn by using the mutations data frame
# in place of the `simulate_seq()` output

# get the driver mutations
driver_mutations <- seq_results$parameters$driver_mutations

# filter germinal mutations
f_seq <- seq_results$mutations %>%
   dplyr::filter(nature!="germinal")

# plotting the VAF histogram filtering out VAFs below 0.02
plot_VAF_histogram(f_seq, cuts = c(0.02, 1))
#> Warning: Missing driver mutations. Disabling driver mutation labelling.


# plotting the VAF histogram with labels
plot_VAF_histogram(f_seq, labels = f_seq["cause"], cuts = c(0.02, 1))
#> Warning: Missing driver mutations. Disabling driver mutation labelling.


# use the driver codes in the driver mutation labels
plot_VAF_histogram(f_seq, labels = f_seq["cause"], cuts = c(0.02, 1),
                   driver_mutations = driver_mutations)


# avoid the driver mutation labels
plot_VAF_histogram(f_seq, labels = f_seq["cause"], cuts = c(0.02, 1),
                   driver_mutation_labels = FALSE)

```
