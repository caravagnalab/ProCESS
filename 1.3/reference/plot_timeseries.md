# Plot the number of cells along a simulation

A pie chart with population counts, split by mutants and epistate. It
also provides annotations for the simulation data.

## Usage

``` r
plot_timeseries(simulation, color_map = NULL)
```

## Arguments

- simulation:

  A simulation.

- color_map:

  A named vector representing the simulation species color map
  (optional).

## Value

A ggplot plot.

## Examples

``` r
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
                                   E2 = 0.01),
                         E2 = list(duplication = 0.08, death = 0.01,
                                   E1 = 0.02)))
sim$history_delta <- 1
sim$place_cell("A[E1]", 500, 500)
sim$run_up_to_time(60)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


plot_timeseries(sim)


# define a custom color map
color_map <- c("#B2DF8A", "#E31A1C")
names(color_map) <- c("A[E1]", "A[E2]")

plot_timeseries(sim, color_map = color_map)
```
