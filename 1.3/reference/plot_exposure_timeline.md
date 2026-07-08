# Plot the signature exposure timeline of a phylogenetic forest

Plots the signatures exposure changes along a phylogenetic forest.

## Usage

``` r
plot_exposure_timeline(
  phylogenetic_forest,
  linewidth = 0.8,
  emphasize_switches = FALSE,
  pal_name = "Dark2"
)
```

## Arguments

- phylogenetic_forest:

  A phylogenetic forest.

- linewidth:

  The width of the lines in the plot.

- emphasize_switches:

  A Boolean flag to emphasize the exposure switches.

- pal_name:

  The name of a `RColorBrewer` palette.

## Value

A `ggplot2` plot.

## Examples

``` r
sim <- TissueSimulation()

sim$add_mutant("A", c(duplication = 0.2))
sim$place_cell("A", 500, 500)

sim$run_up_to_time(150)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# sampling tissue

n_w <- n_h <- 50
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)

sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)

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

 [█████████████████-----------------------] 40% [00m:01s] Processing chr. 22                                                                                                          

 [██████████████████████████████----------] 73% [00m:02s] Processing chr. 22                                                                                                          

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

 [█████████████████████████---------------] 60% [00m:02s] Saving RS index                                                                                                             
done
#> 
 [████████████████████████████████████████] 100% [00m:02s] RS index saved                                                                                                             

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Saving germline                                                                                                              

 [████████████████████████████████████████] 100% [00m:00s] Germline saved                                                                                                             


m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = 8e-8))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                                                                                          

 [█---------------------------------------] 0% [00m:00s] Found 22                                                                                                                     

 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                                                                   

 [███████████-----------------------------] 26% [00m:01s] Reading 22                                                                                                                  

 [████████████████------------------------] 39% [00m:02s] Reading 22                                                                                                                  

 [███████████████████████-----------------] 56% [00m:03s] Reading 22                                                                                                                  

 [██████████████████████████████████------] 84% [00m:04s] Reading 22                                                                                                                  

 [████████████████████████████████████████] 100% [00m:04s] "A"'s SIDs validated                                                                                                       


m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8, ID3 = 1))
m_engine$add_exposure(time = 50,
                      c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5,
                        ID2 = 0.8, ID21 = 0.2))
phylo_forest <- m_engine$place_mutations(forest, 500, 10)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                                                                                            

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                                                                                           


# plotting the phylogenetic forest
plot_exposure_timeline(phylo_forest)


# plotting the phylogenetic forest emphatizing the exposure switches
plot_exposure_timeline(phylo_forest, emphasize_switches = TRUE)


# deleting the mutation engine directory
unlink("demo", recursive = TRUE)
```
