# Plot a tissue

Plots cells distribution over a tissue highlighting species by color. To
facilitate the plot and avoid excessive number of cells, for instance,
when a simulation deals with millions of cells, the plot draws a
hexagonal heatmap of 2D bins.

## Usage

``` r
plot_tissue(
  simulation,
  num_of_bins = 100,
  before_sample = NULL,
  at_sample = NULL,
  plot_next_sample_regions = FALSE,
  plot_sample_region = TRUE,
  color_map = NULL,
  list_all_species = FALSE
)
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
  appeared when the specified sampling occurred. The parameters
  `before_sample` and `at_sample` are mutually exclusive (optional).

- plot_next_sample_regions:

  A Boolean value. When `before_sample` is set and
  `plot_next_sample_regions` is set to be TRUE, this function plots the
  regions of the samples collected at the same simulated time of the
  specified sample. When, instead, `at_sample` is set and
  `plot_next_sample_regions` is set to be TRUE, the function plots the
  regions of the samples collected at the same simulated time of the
  specified sample, but not before the specified sample (default:
  FALSE).

- plot_sample_region:

  A Boolean value. When either `at_sample` or `before_sample` are set
  and `plot_sample_region` is set to be TRUE, the function also plots
  the region of the specified sample (default: TRUE).

- color_map:

  A named vector representing the simulation species color map
  (optional).

- list_all_species:

  A Boolean flag to show all species in the legend (default: FALSE).

## Value

An editable ggplot plot.

## Examples

``` r
set.seed(0)
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
                                   E2 = 0.01),
                         E2 = list(duplication = 0.08, death = 0.01,
                                   E1 = 0.01)))
sim$place_cell("A[E1]", 500, 500)
sim$run_up_to_size("A[E2]", 60000)
#> 
 [█████-----------------------------------] 12% [00m:00s] Cells: 12219                                                                                                           

 [██████████------------------------------] 23% [00m:01s] Cells: 20441                                                                                                           

 [███████████████-------------------------] 37% [00m:02s] Cells: 30792                                                                                                           

 [██████████████████████------------------] 54% [00m:03s] Cells: 43129                                                                                                           

 [████████████████████████████------------] 69% [00m:04s] Cells: 53207                                                                                                           

 [██████████████████████████████████------] 83% [00m:05s] Cells: 62145                                                                                                           

 [███████████████████████████████████████-] 96% [00m:06s] Cells: 71107                                                                                                           

 [████████████████████████████████████████] 100% [00m:06s] Saving snapshot                                                                                                       


# collect 3 samples: "Sample_A", "Sample_B", and "Sample_C"
sim$sample_cells("Sample_A", c(425, 425), c(475, 475))
sim$sample_cells("Sample_B", c(525, 525), c(575, 575))
sim$sample_cells("Sample_C", c(425, 525), c(475, 575))

# let the simulation evolve until the species "A[E2]" account
# for 80000 cells
sim$run_up_to_size("A[E2]", 80000)
#> 
 [███████████████████████████████---------] 75% [00m:00s] Cells: 74726                                                                                                           

 [██████████████████████████████████------] 84% [00m:00s] Cells: 83642                                                                                                           

 [██████████████████████████████████████--] 93% [00m:01s] Cells: 92570                                                                                                           

 [████████████████████████████████████████] 100% [00m:02s] Saving snapshot                                                                                                       


# plot the tissue in the current status
plot_tissue(sim)


# plot the tissue as it was when "Sample_B" was collected
plot_tissue(sim, at_sample = "Sample_B")


# plot the tissue as it was when "Sample_B" was collected and
# highlight the regions of the samples collected at the same
# simulated time, but not before it, i.e., "Sample_B" and
# "Sample_C"
plot_tissue(sim, at_sample="Sample_B",
            plot_next_sample_regions = TRUE)


# plot the tissue as it was just before sampling "Sample_B"
plot_tissue(sim, before_sample="Sample_B")


# plot the tissue as it was just before sampling "Sample_B"
# and highlight the regions of the samples collected at the
# same simulated time, i.e., "Sample_A", "Sample_B", and
# "Sample_C"
plot_tissue(sim, before_sample="Sample_B",
            plot_next_sample_regions = TRUE)


# define a custom color map
color_map <- c("#B2DF8A", "#E31A1C")
names(color_map) <- c("A[E1]", "A[E2]")

plot_tissue(sim, color_map = color_map)
```
