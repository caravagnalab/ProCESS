# Getting forest species

This method builds a data frame containing information about the
simulated species.

## Value

A data frame reporting `mutant` and, if the simulation has epigenetic
states, `epistate` for each registered species.

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add the mutant "A"
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation until "A" has less than 15 cells
sim$run_up_to_size("A", 15)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# add the mutant "B"
sim$add_mutant("B", c(duplication = 0.3, death = 0.01))

# let one border cell of "A" generate a cell in "B"
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# run the simulation until "B" has less than 100 cells
sim$run_up_to_size("B", 30)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# add the mutant "C"
sim$add_mutant("C", c(duplication = 0.4, death = 0.01))

# let one border cell of "B" generate a cell in "C"
sim$mutate_progeny(sim$choose_border_cell_in("B"), "C")

# run the simulation until "C" has less than 2000 cells
sim$run_up_to_size("C", 2000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# search for a 33x33 region containing 50 cells in A and
# 50 cells in B at least and sample it
region <- sim$search_sample(c(A = 50, B = 50), 33, 33)
sim$sample_cells("S1", region$lower_corner, region$upper_corner)

# search for a 33x33 region containing 50 cells in B and
# 50 cells in C at least and sample it
region <- sim$search_sample(c(B = 50, C = 50), 33, 33)
sim$sample_cells("S2", region$lower_corner, region$upper_corner)

# build the sample forest
forest <- sim$get_sample_forest()

# get species information. Since the simulation has no epigenetic
# state, the species correspond to the mutants
forest$get_species_info()
#>   mutant
#> 1      A
#> 2      B
#> 3      C

# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add the mutant "A"
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.01,
                                   E2 = 0.05),
                         E2 = list(duplication = 0.2, death = 0.01,
                                   E2 = 0.01)))

# place a cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# run the simulation until "A" has less than 15 cells
sim$run_up_to_size("A[E2]", 40)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# add the mutant "B"
sim$add_mutant("B", list(E2 = list(duplication = 0.3,
                                   death = 0.01)))

# let one border cell of "A[E2]" generate a cell in "B"
sim$mutate_progeny(sim$choose_border_cell_in("A[E2]"), "B")

# run the simulation until "B[E2]" has less than 100 cells
sim$run_up_to_size("B[E2]", 50)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# add the mutant "C"
sim$add_mutant("C", list(E2 = list(duplication = 0.5,
                                   death = 0.01)))

# let one border cell of "B" generate a cell in "C"
sim$mutate_progeny(sim$choose_border_cell_in("B"), "C")

# run the simulation until "C" has less than 2000 cells
sim$run_up_to_size("C[E2]", 2000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# search for a 33x33 region containing 50 cells in A and
# 50 cells in B at least and sample it
region <- sim$search_sample(c(A = 50, B = 50), 33, 33)
sim$sample_cells("S1", region$lower_corner, region$upper_corner)

# search for a 33x33 region containing 50 cells in B and
# 50 cells in C at least and sample it
region <- sim$search_sample(c(B = 50, C = 50), 33, 33)
sim$sample_cells("S2", region$lower_corner, region$upper_corner)

# build the sample forest
forest <- sim$get_sample_forest()

# get species information
forest$get_species_info()
#>   mutant epistate
#> 1      A       E1
#> 2      A       E2
#> 3      B       E1
#> 4      B       E2
#> 5      C       E1
#> 6      C       E2
```
