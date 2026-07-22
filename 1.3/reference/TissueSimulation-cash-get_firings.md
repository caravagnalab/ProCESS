# Getting the number of fired events

This method returns a data frame reporting the current number of
simulated events per species.

## Value

A data frame having the `event`, `mutant`, `fired`, and, when the
simulation has epigenetic states, `epistate`. Each row reports event of
a given type have been fired in a species.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.1))

# add mutant "B" and set its rates
sim$add_mutant("B", c(duplication = 0.15, death = 0.05))

# schedule an evolution from mutant "A" to mutant "B" at time 10
sim$schedule_mutation("A", "B", 10)

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                       


# get the history of species counts
sim$get_firings()
#>          event mutant fired
#> 1       deaths      A   936
#> 2 duplications      A  1523
#> 3       deaths      B     8
#> 4 duplications      B    11
# set the seed of the random number generator
set.seed(0)

# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# add mutant "B" and set its species rates
sim$add_mutant("B",
               list(E1 = list(duplication = 0.4, death = 0.1, E2 = 0.02),
                    E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))

# schedule a mutation from "A" to "B"
sim$schedule_mutation("A", "B", 10)

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                       


# counts the number of cells per species
sim$get_firings()
#>           event mutant epistate fired
#> 1        deaths      A       E1   555
#> 2  duplications      A       E1   719
#> 3      switches      A       E1    31
#> 4        deaths      A       E2    20
#> 5  duplications      A       E2    87
#> 6      switches      A       E2     6
#> 7        deaths      B       E1  3538
#> 8  duplications      B       E1  6484
#> 9      switches      B       E1   307
#> 10       deaths      B       E2    52
#> 11 duplications      B       E2   218
#> 12     switches      B       E2    21
```
