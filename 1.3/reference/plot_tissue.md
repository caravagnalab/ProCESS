# Plot a tissue

This function plots the tissue

## Usage

``` r
plot_tissue(simulation)
```

## Arguments

- simulation:

  A simulation object.

- num_of_bins:

  The number of bins (default: 100).

- before_sample:

  A sample name. When provided, this function represents the tissue as
  appeared just before the specified sampling. The parameters
  `before_sample` and `at_sample` are mutually exclusive (optional).

- at_sample:

  A sample name. When provided, this function represents the tissue as
  appeared when the specified sampling was about to be collected. The
  parameters `before_sample` and `at_sample` are mutually exclusive
  (optional).

- plot_next_sample_regions:

  A Boolean value. When `before_sample` is set and
  `plot_next_sample_regions` is set to be `TRUE`, this function plots
  the regions of the samples collected at the same simulated time of the
  specified sample. When, instead, `at_sample` is set and
  `plot_next_sample_regions` is set to be `TRUE`, the function plots the
  regions of the samples collected at the same simulated time of the
  specified sample, but not before the specified sample (default:
  `FALSE`).

- plot_sample_region:

  A Boolean value. When either `at_sample` or `before_sample` are set
  and `plot_sample_region` is set to be `TRUE`, the function also plots
  the region of the specified sample (default: `TRUE`).

- color_map:

  A named vector representing the simulation species color map
  (optional).

- list_all_species:

  A Boolean flag to show all species in the legend (default: `FALSE`).

- focus_function:

  A function whose input is the data frame returned by
  [`TissueSimulation$get_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_cells.md)
  and returns a Boolean vector whose length is the number of rows in the
  input data frame. When one the row in the output is `FALSE` the
  corresponding cells is plotted in grey. When the parameter is set to
  `NULL`, all tumour simulation cells are colored (default: `NULL`).

- alpha_function:

  A function whose input is the data frame returned by
  [`TissueSimulation$get_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_cells.md)
  and returns a real vector whose values are in the interval \\\[0,1\]\\
  and whose length is the number of rows in the input data frame. Each
  value in the output is used as alpha level of the corresponding cell.
  When the parameter is set to `NULL`, all tumour simulation cells have
  alpha level `1` (default: `NULL`).

## Value

An editable ggplot plot.

## Details

This function plots cells distribution over a tissue highlighting
species by color. To facilitate the plot and avoid excessive number of
cells, for instance, when a simulation deals with millions of cells, the
plot draws a hexagonal heatmap of 2D bins.

## See also

[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)

## Examples

``` r
# set the seed
set.seed(0)

# build a tissue simulation
sim <- TissueSimulation(width = 600, height = 600)

# avoid drift
sim$death_activation_level <- 50

# add the mutant A
sim$add_mutant("A", c(duplication = 0.12, death = 0.05))

# place a cell in the tissue and simulate it until 10 cells
sim$place_cell("A", 300, 300)
sim$run_up_to_size("A", 10, quiet = TRUE)

# add the mutant B and let mutate a border cell of A in B
sim$add_mutant("B", c(duplication = 0.145, death = 0.06))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# simulate the tissue up to 30 cells in B
sim$run_up_to_size("B", 30, quiet = TRUE)

# add the third mutant and let one cell of A mutate into C
sim$add_mutant("C", c(duplication = 0.15, death = 0.06))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "C")

# simulate the tissue until C consists of 25000 cells
sim$run_up_to_size("C", 25000, quiet = TRUE)

# collect the sample "S1"
sim$sample_cells("S1", c(145, 230), c(215, 300))

# let the simulation reach 25000 cells in C again
sim$run_up_to_size("C", 25000, quiet = TRUE)

# collect two samples
sim$sample_cells("S2", c(350, 300), c(420, 370))
sim$sample_cells("S3", c(200, 350), c(270, 420))

# add a further mutant and derive it from B
sim$add_mutant(name = "D", c(duplication = 0.8, death = 0.01))
sim$mutate_progeny(sim$choose_border_cell_in("B"), "D")

# let the tumour evolve until the mutant C and D cumulatively
# consist of 10000 cells
sim$run_until(sim$var("C") + sim$var("D") == 1e5, 10, quiet = TRUE)

# plot the tissue in the current status
plot_tissue(sim)


# plot the tissue as it was when "S3" was about to be sampled
plot_tissue(sim, at_sample = "S3")


# plot the tissue as it was when "S3" was about to be sampled and
# highlight the regions of the samples collected at the same simulated
# time, but not before
# it, i.e., "S3"
plot_tissue(sim, at_sample="S3",
            plot_next_sample_regions = TRUE)


# plot the tissue as it was just before sampling "S3"
plot_tissue(sim, before_sample="S3")


# plot the tissue as it was just before sampling "S3" and highlight
# the regions of the samples collected at the same simulated time,
# i.e., "S2" and "S3"
plot_tissue(sim, before_sample="S3",
            plot_next_sample_regions = TRUE)


# define a custom color map
color_map <- c(A="#B2DF8A", B="#E31A1C", C="#C41E4E", D="#FEAAAA")
names(color_map) <- c("A", "B", "C", "D")

plot_tissue(sim, color_map = color_map)


# this function returns `TRUE` for cells in the rectangle
# [250,350]x[300,350]
focus_function <- function(cells) {
  (cells$position_x >= 250 & cells$position_x <= 350
   & cells$position_y >= 300 & cells$position_y <= 350)
}

# plot the tissue highlighting the region [250,350]x[300,350]
plot_tissue(sim, focus_function = focus_function)
```
