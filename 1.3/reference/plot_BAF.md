# Plot the genome-wide B-Allele Frequency (BAF)

This function plots the genome-wide B-Allele Frequency (BAF) of a
specific sample.

## Usage

``` r
plot_BAF(
  seq_result,
  chromosomes = NULL,
  samples = NULL,
  labels = NULL,
  mutation_filter = filter_germinal,
  driver_mutations = NULL,
  driver_mutation_labels = TRUE,
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
  `function(x) x %>% dplyr::filter(classes != "germinal")`).

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

- N:

  The number of mutations to sample for plotting (default: 5000).

## Value

A ggplot2 object showing the BAF distribution across the genome.

## See also

[`plot_VAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF.md),
[`plot_DR()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_DR.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.2))
sim$place_cell("A", 500, 500)
sim$run_up_to_time(50)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sampling tissue
n_w <- n_h <- 10
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("Sample_A", bbox$lower_corner, bbox$upper_corner)

# adding second mutant
sim$add_mutant("B", c(duplication = 0.3))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
sim$run_up_to_time(300)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sampling tissue again
bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
sim$sample_cells("Sample_B", bbox$lower_corner, bbox$upper_corner)

forest <- sim$get_sample_forest()

# placing mutations
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                     

 [█████████████---------------------------] 30% [00m:01s] Loading RS index                                    

 [█████████████████████████---------------] 61% [00m:02s] Loading RS index                                    

 [████████████████████████████████████----] 89% [00m:03s] Loading RS index                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                    


m_engine$add_mutant(mutant_name="A", passenger_rates = c(SNV = 5e-8),
                    drivers = list(SNV("22", 46510210, "C", "A", allele = 1),
                                   "DGCR8 P26L"))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                  

 [████████████████████████████████████████] 100% [00m:00s] "A"'s SIDs validated                               

m_engine$add_mutant(mutant_name="B", passenger_rates = c(SNV = 5e-9),
                    drivers = list(list("DGCR8 A18V", allele = 1)))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs                                  

 [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated                               

m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))

phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                    

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                   


# simulate sequencing and avoid progress bar
seq_results <- simulate_seq(phylo_forest, coverage = 10,
                            write_SAM = FALSE, quiet = TRUE)

# plotting the BAF over all samples
plot_BAF(seq_results)


# plotting the BAF for Sample_B only
plot_BAF(seq_results, samples = c("Sample_B"))


# plotting the BAF for Sample_B with labels
plot_BAF(seq_results, samples = "Sample_B",
         labels = seq_results$mutations["classes"])


# let us define a function to filter germinal and pre-neoplastic
# from the input data
library(dplyr)
filter_data <- function(data) {
  data %>% dplyr::filter(!classes %in% list("germinal",
                                            "pre-neoplastic"))
}

# plotting the BAF without germinal and pre-neoplastic
plot_BAF(seq_results, mutation_filter = filter_data)


# the same plots can be drawn by using the mutations data frame
# in place of the `simulate_seq()` output

# filter germinal mutations
f_seq <- seq_results$mutations %>%
   dplyr::filter(classes != "germinal")

# plotting the BAF
plot_BAF(f_seq)
#> Warning: Missing driver mutations. Disabling driver mutation labelling.


unlink('demo', recursive = T)
```
