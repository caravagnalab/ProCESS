# Retrieving sample information

This method retrieves information about the samples collected along the
simulation.

## Value

A data frame containing, for each sample collected during the
simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`, `ymax`,
`xmax`, `tumour_cells`, and `tumour_cells_in_bbox`. The columns `ymin`,
`xmin`, `ymax`, `xmax` report the boundaries of the sample bounding box,
while `tumour_cells` and `tumour_cells_in_bbox` are the number of tumour
cells in the sample and in the bounding box, respectively.

## See also

[`TissueSimulation$sample_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-sample_cells.md),
[`SampleForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_samples_info.md),
[`PhylogeneticForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_samples_info.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", list(duplication = 0.3, death = 0.01))

# place an "A" cell in the tissue
sim$place_cell("A", 500, 500)

# simulate the tissue until "A" consists of 50000 cells
sim$run_up_to_size("A", 50000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 475),
                 upper_corner = c(500, 550))

# simulate 1 time unit more
sim$run_up_to_time(sim$get_clock()+1)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sample the region [500,520]x[525,550]
sim$sample_cells("S2", lower_corner = c(500, 525),
                 upper_corner = c(520, 550))

# get information about all the collected
# samples, i.e, S1 and S2
sim$get_samples_info()
#>   name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1   S1 17  450  475  500  550         3851                 3851 197.9530
#> 2   S2 18  500  525  520  550          520                  520 198.9531
```
