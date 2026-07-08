# Searching rectangular tissue samples

This method searches a set of rectangular tissue samples.

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

- n_samples:

  The number of searched samples.

- seed:

  The seed of the random generator the select the samples among those
  satisfying the constraints (optional).

## Value

A vector of `n_samples` rectangular tissue samples that satisfy the
aimed constraints.

## Details

The aimed samples mush satisfy the specified number of cells. The sizes
of the samples are also provided a parameter of the method. This method
takes asymptotic time \\O(\|\textrm{tissue width}\|\*\|\textrm{tissue
height}\|)\\.

## See also

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation(width = 150, height = 150)

# add mutant "A" and set its rates
sim$add_mutant("A", list(duplication = 0.20, death = 0.01))

# place an "A" cell in the tissue
sim$place_cell("A", 75, 75)

# simulate the tissue until "A" consists of 300 cells
sim$run_up_to_size("A", 130)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# add mutant "B" and set its rates
sim$add_mutant("B", list(duplication = 0.3, death = 0.01))

# mutate a border cell in "A" into "B"
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# run the simulation until "B" consists of 1000 cells
sim$run_up_to_size("B", 4000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# plot the tissue as it is
plot <- plot_tissue(sim)

# find 3 50x50 samples containing 80 "B" cells and 100 "A" cells
# at least
bboxes <- sim$search_samples(min_num_of_cells = c("A" = 100, "B" = 80),
                             width = 25, height = 25,
                             n_samples = 3)
bboxes
#> [[1]]
#> TissueRectangle((31,63),(55,87))
#> 
#> [[2]]
#> TissueRectangle((31,88),(55,112))
#> 
#> [[3]]
#> TissueRectangle((56,63),(80,87))
#> 

# plot the found bounding boxes
for (bbox in bboxes) {
  plot <- plot +
    ggplot2::geom_rect(xmin = bbox$lower_corner[1],
                       xmax = bbox$upper_corner[1],
                       ymin = bbox$lower_corner[2],
                       ymax = bbox$upper_corner[2],
                       fill = NA, color = "black")
}

plot
```
