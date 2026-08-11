# Searching for a rectangular tissue sample

This method searches a rectangular tissue sample.

## Arguments

- min_num_of_cells:

  A named integer vector reporting the minimum number of cells per
  species or mutant.

- num_of_cells:

  The number of cells in the searched sample.

- width:

  The width of the searched sample.

- height:

  The height of the searched sample.

## Value

If a rectangular sample satisfying the provided constraints can be
found, the corresponding rectangle.

## Details

The aimed sample must satisfy the specified number of cells. The size of
the sample is also provided a parameter of the method. The complexity of
this method is \\O(\|\textrm{tissue width}\|\*\|\textrm{tissue
height}\|)\\.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

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

# simulate the tissue until "A" consists of 100 cells
sim$run_up_to_size("A", 100)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# add mutant "B" and set its rates
sim$add_mutant("B", list(duplication = 0.3, death = 0.01))

# mutate a border cell in "A" into "B"
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# run the simulation until "B" consists of 1000 cells
sim$run_up_to_size("B", 1000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# find a 50x50 sample containing 80 "B" cells and 10 "A" cells at least
sim$search_sample(c("A" = 10, "B" = 80), 50, 50)
#> TissueRectangle((462,454),(511,503))
```
