# Plot Marginals of Variant Allele Frequency (VAF)

This function generates scatter plots showing the marginal distributions
of Variant Allele Frequency (VAF) for pairs of samples on a specific
chromosome.

## Usage

``` r
plot_VAF_marginals(
  seq_result,
  chromosomes = NULL,
  samples = NULL,
  labels = NULL,
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
  plot. A mutation is represented in the plot only if its VAF lays in
  this interval in at least one of the samples (default: `c(0, 1)`).

## Value

A list of ggplot2 objects showing scatter plots of VAF marginals for
pairs of samples.

## See also

[`plot_VAF_histogram()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF_histogram.md),
[`plot_VAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF.md)

## Examples

``` r
# use a sequencing result example
seq_results <- example("Sequencing results")

# plotting the VAF marginals without germinal mutations
plot_VAF_marginals(seq_results)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# let us define a function to filter germinal and pre-neoplastic
# from the input data
library(dplyr)
filter_data <- function(data) {
  data %>% dplyr::filter(!nature %in% list("germinal",
                                           "pre-neoplastic"))
}

# plotting the VAF marginals without germinal and pre-neoplastic
plot_VAF_marginals(seq_results, mutation_filter = filter_data)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# plotting the VAF marginals filtering out mutations having
# the VAF below 0.2 in both the samples
plot_VAF_marginals(seq_results, cuts = c(0.2, 1))
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# plotting the VAF marginals with labels
plot_VAF_marginals(seq_results, labels = seq_results$mutations["cause"])
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# avoid the driver mutation labels
plot_VAF_marginals(seq_results, labels = seq_results$mutations["cause"],
                   driver_mutation_labels = FALSE)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# the same plots can be drawn by using the mutations data frame
# in place of the `simulate_seq()` output

# get driver mutations
driver_mutations <- seq_results$parameters$driver_mutations

# filter germinal mutations
f_seq <- seq_results$mutations %>%
   dplyr::filter(nature!="germinal")

# plotting the VAF histogram filtering out VAFs below 0.02
plot_VAF_marginals(f_seq, cuts = c(0.2, 1))
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# plotting the VAF histogram with labels
plot_VAF_marginals(f_seq, labels = f_seq["cause"])
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# use the driver codes in the driver mutation labels
plot_VAF_marginals(f_seq, labels = f_seq["cause"],
                   driver_mutations = driver_mutations)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 

# avoid the driver mutation labels
plot_VAF_marginals(f_seq, labels = f_seq["cause"],
                   driver_mutation_labels = FALSE)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
#> [[4]]

#> 
#> [[5]]

#> 
#> [[6]]

#> 
```
