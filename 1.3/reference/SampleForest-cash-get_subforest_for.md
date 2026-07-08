# Building sub-forests

This method builds a sub-forest using as leaves some of the original
samples.

## Arguments

- sample_names:

  The names of the samples whose cells will be used as leaves of the new
  forest.

## Value

A sample forest built on the samples mentioned in `sample_names`

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add the mutant "A"
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation until "A" has less than 50000 cells
sim$run_up_to_size("A", 50000)
#> 
 [██████████████████████████████████------] 83% [00m:00s] Cells: 41559                                        

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 500), upper_corner = c(475, 550))

# sample the region [550,650]x[600,675]
sim$sample_cells("S2", lower_corner = c(550, 600), upper_corner = c(650, 675))

# build the sample forest
forest <- sim$get_sample_forest()

# show the forest data
forest
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 5765
#>   # of leaves: 1561
#>   samples: {"S1", "S2"}
#> 

# get the subforest for sample "S2"
forest$get_subforest_for("S2")
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 868
#>   # of leaves: 240
#>   samples: {"S2"}
#> 
```
