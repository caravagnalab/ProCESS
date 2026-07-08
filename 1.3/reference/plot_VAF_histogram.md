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

## Value

A ggplot2 object showing the VAF histogram.

## See also

[`plot_VAF_marginals()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF_marginals.md),
[`plot_VAF()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.1))
sim$place_cell("A", 500, 500)
sim$run_up_to_time(100)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# sampling tissue
n_w <- n_h <- 10
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("SampleA", bbox$lower_corner, bbox$upper_corner)

# adding second mutant
sim$add_mutant("B", c(duplication = 0.3))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
sim$run_up_to_time(300)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# sampling tissue again
bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
sim$sample_cells("SampleB", bbox$lower_corner, bbox$upper_corner)

forest <- sim$get_sample_forest()

# placing mutations
m_engine <- MutationEngine(setup_code = "demo")
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
#> Downloading signature files...
#> Signature file downloaded
#> Downloading driver mutation file...
#> Driver mutation file downloaded
#> Decompressing driver mutation file...done
#> Downloading passenger CNAs file...
#> Passenger CNAs file downloaded
#> Decompressing passenger CNAs file...done
#> Downloading germline...
#> Germline downloaded
#> Decompressing mutations...
#> done
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                                                           

 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22                                                                                                          

 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22                                                                                                          

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                                                        

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                                                         

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                                                        

#> done
#> Building repeated sequence index...
#> 
 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                                                                   

 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                                                                 

#> 
 [████████████████████████████████████████] 100% [00m:01s] RS index built                                                                                                             

#> 
 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                                                              

 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                                                              

 [████████████████------------------------] 38% [00m:02s] Saving RS index                                                                                                             
done
#> 
 [████████████████████████████████████████] 100% [00m:02s] RS index saved                                                                                                             

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Saving germline                                                                                                              

 [████████████████████████████████████████] 100% [00m:00s] Germline saved                                                                                                             


m_engine$add_mutant("A", passenger_rates = c(SNV = 5e-8),
                    drivers = list(SNV("22", 16510210, "C", "T", allele = 1),
                                   "DGCR8 P26L"))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                                                                                          

 [█---------------------------------------] 0% [00m:00s] Found 22                                                                                                                     

 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                                                                   

 [█████████████---------------------------] 30% [00m:01s] Reading 22                                                                                                                  

 [█████████████████████████---------------] 60% [00m:02s] Reading 22                                                                                                                  

 [█████████████████████████████████████---] 92% [00m:03s] Reading 22                                                                                                                  

 [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated                                                                                                       

m_engine$add_mutant(mutant_name = "B", passenger_rates = c(SNV = 5e-9),
                    drivers = list("DGCR8 A18V"))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs                                                                                                          

 [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated                                                                                                       

m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))

phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                                                                                            

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                                                                                           


# simulating sequencing without the normal sample
seq_results <- simulate_seq(phylo_forest, coverage = 10, write_SAM = F,
                            with_normal_sample = FALSE, quiet = TRUE)

# plotting the VAF histogram without germinal mutations
plot_VAF_histogram(seq_results)


# let us define a function to filter germinal and pre-neoplastic
# from the input data
library(dplyr)
filter_data <- function(data) {
  data %>% dplyr::filter(!classes %in% list("germinal",
                                            "pre-neoplastic"))
}

# plotting the VAF histogram without germinal and pre-neoplastic
plot_VAF_histogram(seq_results, mutation_filter = filter_data)


# plotting the VAF histogram filtering out VAFs below 0.02
plot_VAF_histogram(seq_results, cuts = c(0.02, 1))


# plotting the VAF histogram with labels
plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
                   labels = seq_results$mutations["causes"])


# avoid the driver mutation labels
plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
                   labels = seq_results$mutations["causes"],
                   driver_mutation_labels = FALSE)


# the same plots can be drawn by using the mutations data frame
# in place of the `simulate_seq()` output

# filter germinal mutations
f_seq <- seq_results$mutations %>%
   dplyr::filter(classes!="germinal")

# plotting the VAF histogram filtering out VAFs below 0.02
plot_VAF_histogram(f_seq, cuts = c(0.02, 1))
#> Warning: Missing driver mutations. Disabling driver mutation labelling.


# plotting the VAF histogram with labels
plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1))
#> Warning: Missing driver mutations. Disabling driver mutation labelling.


# use the driver codes in the driver mutation labels
plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1),
                   driver_mutations = phylo_forest$get_driver_mutations())


# avoid the driver mutation labels
plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1),
                   driver_mutation_labels = FALSE)


# deleting the mutation engine directory
unlink('demo', recursive = T)
```
