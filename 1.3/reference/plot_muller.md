# Draw a Muller plot of the simulation

Plots a muller plot of the simulation, using the store time-series data
and package `ggmuller`.

## Usage

``` r
plot_muller(simulation, color_map = NULL)
```

## Arguments

- simulation:

  A simulation object.

- color_map:

  A named vector representing the simulation species color map
  (optional).

## Value

A `ggplot2` plot.

## Examples

``` r
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.02,
                                   E2 = 0.01),
                         E2 = list(duplication = 0.08, death = 0.01,
                                   E1 = 0.01)))
sim$history_delta <- 1
sim$place_cell("A[E1]", 500, 500)
sim$run_up_to_time(60)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot             


plot_muller(sim)


# define a custom color map
color_map <- c("#B2DF8A", "#E31A1C")
names(color_map) <- c("A[E1]", "A[E2]")

plot_muller(sim, color_map = color_map)
```
