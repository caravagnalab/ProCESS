# Getting the sample forest

This method returns the sample forest.

## Value

The sample forest having as leaves the sampled cells

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
sim$place_cell("A", 500, 500)

sim$death_activation_level <- 100
sim$run_up_to_size("A", 50000)
#> 
 [███████████████████████████████████████-] 96% [00m:00s] Cells: 48036                                                                        

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                    


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", c(450,475), c(500,550))

# build the sample forest
forest <- sim$get_sample_forest()

forest
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 13489
#>   # of leaves: 3849
#>   samples: {"S1"}
#> 
```
