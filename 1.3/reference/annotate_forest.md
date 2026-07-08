# Annotate a plot of cell divisions

It annotates a plot of cell divisions with information from sampling
times and MRCAs for all available samples

## Usage

``` r
annotate_forest(
  tree_plot,
  forest,
  samples = TRUE,
  MRCAs = TRUE,
  exposures = FALSE,
  facet_signatures = TRUE,
  drivers = TRUE,
  add_driver_label = TRUE
)
```

## Arguments

- tree_plot:

  The output of `plot_forest`.

- forest:

  The original forest object from which the input to `plot_forest`has
  been derived.

- samples:

  If `TRUE` it annotates samples.

- MRCAs:

  If `TRUE` it annotates MRCAs.

- exposures:

  If `TRUE` it annotates exposures to mutational signatures.

- facet_signatures:

  If `TRUE` and if `exposures` is `TRUE` it creates a faceted forest
  plot where the exposure to each signature is annotated on a separated
  plot.

- drivers:

  If `TRUE` it annotates drivers on the node they originated.

- add_driver_label:

  If `TRUE` and if `drivers` is `TRUE` it annotates the driver name.

## Value

A `ggraph` tree plot.

## Examples

``` r
sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.08, death = 0.01))
sim$place_cell("A", 500, 500)
sim$run_up_to_time(60)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    

sim$sample_cells("MySample", c(500, 500), c(510, 510))
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


m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = 1e-9),
                    drivers = list(SNV("22", 10510210, "C"),
                                   CNA(type = "A", "22", from = 10303470,
                                       len = 200000)))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                  

 [████████████████████████████████████████] 100% [00m:00s] "A"'s SIDs validated                               

m_engine$add_exposure(c(SBS13 = 0.2, SBS1 = 0.8))
m_engine$add_exposure(time = 50, c(SBS17b = 0.2, SBS3 = 0.8))

forest <- sim$get_sample_forest()
forest$get_samples_info()
#>       name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 MySample 22  500  500  510  510           27                   27 60.39314
forest_muts <- m_engine$place_mutations(forest, 1000, 500)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                    

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                   

tree_plot <- plot_forest(forest)
annotate_forest(tree_plot, forest_muts, samples = T, MRCAs = T,
                exposures = T, drivers = T, add_driver_label = T)
```
