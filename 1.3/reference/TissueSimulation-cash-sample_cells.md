# Sampling a set of cells

This method samples a set of tumour cells.

## Arguments

- sample_name:

  The name of the sample.

- lower_corner:

  The lower corner of the sample bounding box (optional in pair with
  `upper_corner`).

- upper_corner:

  The upper corner of the sample bounding box (optional in pair with
  `lower_corner`).

- num_of_cells:

  The maximum number of tumour cells to collect (optional).

## Details

It removes the cells from the simulated tissue and stores them in a
sample that can be subsequently retrieved to build a sample forest.

## See also

[`TissueSimulation$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", list(duplication = 0.2, death = 0.01))

# place an "A" cell in the tissue
sim$place_cell("A", 500, 500)

# simulate the tissue until "A" consists of 50000 cells
sim$run_up_to_size("A", 50000)
#> 
 [████████████████████████████████████----] 89% [00m:00s] Cells: 44941                                                                                                                

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# randomly sample 50 tumour cells from the tissue
sim$sample_cells(sample_name = "S1", num_of_cells = 50)

# sample the region [450,500]x[475,550]
sim$sample_cells(sample_name = "S2",
                 lower_corner = c(450, 475), upper_corner = c(500, 550))

# randomly sample 50 tumour cells from the tissue region [500,550]x[500,550]
sim$sample_cells(sample_name = "S3",
                 lower_corner = c(500, 500), upper_corner = c(550, 550),
                 num_of_cells = 50)

sim$get_samples_info()
#>   name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1   S1 19    0    0  999  999           50                50000 294.3848
#> 2   S2 20  450  475  500  550         3849                 3849 294.3848
#> 3   S3 21  500  500  550  550           50                 2534 294.3848
```
