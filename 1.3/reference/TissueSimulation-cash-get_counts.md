# Counting the cell number

This method returns the current number of cells per species and that
since the simulation began.

## Value

A data frame reporting `mutant`, `counts`, `overall`, and, when the
simulation has epigenetic states, `epistate`. The data frame contains a
row for each species in the simulation.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# add mutant "B" and set its rates
sim$add_mutant("B", c(duplication = 0.4, death = 0.01))

# schedule an evolution from mutant "A" to mutant "B" at time 20
sim$schedule_mutation("A", "B", 20)

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# counts the number of cells per species
sim$get_counts()
#>   mutant counts overall
#> 1      A    913    2253
#> 2      B   2857    6453
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
sim$get_counts()
#>   mutant epistate counts overall
#> 1      A       E1    139    1445
#> 2      A       E2     92     205
#> 3      B       E1   2661   12990
#> 4      B       E2    452     743
```
