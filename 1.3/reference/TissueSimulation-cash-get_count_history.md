# Getting the history of the number of cells per species

This method returns a data frame reporting the number of species cells
in each sampled simulation time.

## Value

A data frame reporting `mutant`, `counts`, and `time` for each species,
and for each sampled time. When the simulation has epigenetic states,
the data frame also contains the column `epistate`.

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

# set delta time between species counting to 20
sim$history_delta <- 20

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# get the history of species counts
sim$get_count_history()
#>   mutant count     time
#> 1      A     8 20.00000
#> 2      B     1 20.00000
#> 3      A    98 40.00000
#> 4      B     1 40.00000
#> 5      A   359 60.00000
#> 6      B     4 60.00000
#> 7      A   587 70.03108
#> 8      B     4 70.03108
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

# set delta time between species counting to 30
sim$history_delta <- 30

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# counts the number of cells per species
sim$get_count_history()
#>    mutant epistate count    time
#> 1       A       E1    38 30.0000
#> 2       A       E2    11 30.0000
#> 3       B       E1    23 30.0000
#> 4       B       E2     2 30.0000
#> 5       A       E1   174 60.0000
#> 6       A       E2    71 60.0000
#> 7       B       E1  1406 60.0000
#> 8       B       E2   226 60.0000
#> 9       A       E1   139 70.0011
#> 10      A       E2    92 70.0011
#> 11      B       E1  2661 70.0011
#> 12      B       E2   452 70.0011
```
