# Plot the number of stochastic events in the simulation

A pie chart with events split by type, mutant and epigentic state where
they occurred. It also provides annotations for the simulation
information.

## Usage

``` r
plot_firings(simulation)
```

## Arguments

- simulation:

  A simulation.

## Value

A ggplot plot.

## Examples

``` r
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
                                   E2 = 0.01),
                         E2 = list(duplication = 0.08, death = 0.01,
                                   E1 = 0.02)))
sim$place_cell("A[E1]", 500, 500)
sim$run_up_to_time(60)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    

plot_firings(sim)
```
