# Getting the simulation lineage graph

This method returns the lineage graph of the simulation.

## Value

A data frame reporting `ancestor`, `progeny`, and `first_occurrence` of
each species-to-species transition.

## Details

At the beginning of the computation only the species of the added cells
are present in the tissue. As the simulation proceeds new species arise
as a consequence of either mutant mutations or epigenetic switches. The
*lineage graph* stores these species evolutions and it reports the first
occurrence time of any species-to-species transition.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# add mutant "B" and set its species rates
sim$add_mutant("B",
               list(E1 = list(duplication = 0.3, death = 0.1, E2 = 0.02),
                    E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))

# schedule a mutation from "A" to "B"
sim$schedule_mutation("A", "B", 20)

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                    


sim$get_lineage_graph()
#>    ancestor progeny first_cross
#> 1 Wild-type   A[E1]    0.000000
#> 2     A[E1]   A[E2]    5.726872
#> 3     A[E1]   B[E1]   20.000000
#> 4     A[E2]   A[E1]   33.073395
```
