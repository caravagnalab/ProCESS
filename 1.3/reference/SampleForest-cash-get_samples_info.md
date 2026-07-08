# Retrieving the samples' information

This method retrieves information about the samples whose cells were
used as leaves of the sample forest.

## Value

A data frame containing, for each sample collected during the
simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`, `ymax`,
`xmax`, `tumour_cells`, and `tumour_cells_in_bbox`. The columns `ymin`,
`xmin`, `ymax`, `xmax` report the boundaries of the sample bounding box,
while `tumour_cells` and `tumour_cells_in_bbox` are the number of tumour
cells in the sample and in the bounding box, respectively.

## See also

[`PhylogeneticForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_samples_info.md)
for usage examples,
[`TissueSimulation$sample_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-sample_cells.md),
[`TissueSimulation$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)

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
 [█████████████████████████████████-------] 80% [00m:00s] Cells: 40474                                                                                                                

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 500), upper_corner = c(475, 550))

# build the sample forest
forest <- sim$get_sample_forest()

# get information about the sampled whose cells
# are the forest leaves, i.e, S1 and S2
forest$get_samples_info()
#>   name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1   S1  7  450  500  475  550         1321                 1321 294.3848
```
