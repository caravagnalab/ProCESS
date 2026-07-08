# Computing the forest sticks

This method computes the forest sticks.

## Arguments

- birth_threshold:

  The maximum birth time for the cells associated to the returned sticks
  (optional).

## Value

The list of the forest sticks whose associated cells have birth time
smaller than or equal to `birth_threshold`. Each stick is represented as
the list of cell identifiers labelling the nodes in the stick from the
higher to the deeper in the forest.

## Details

A *crucial node* of a forest is a root of the forest, a node whose
parent belongs to a different species, or the most recent common
ancestor of two crucial nodes.

A *stick* is a path of the forest in which the only crucial nodes are
the first and the last one.

This method returns the list of the forest sticks. Each stick is
represented by the sequence of cell identifiers labelling the nodes in
the stick.

## See also

[`PhylogeneticForest$get_sticks()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sticks.md)

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


sim$get_clock()
#> [1] 15.26943

# add the mutant "B"
sim$add_mutant("B", c(duplication = 0.3, death = 0.01))

# let one border cell of "A" generate a cell in "B"
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# run the simulation until "B" has less than 100 cells
sim$run_up_to_size("B", 30)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


sim$get_clock()
#> [1] 49.39972

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

# search for the forest sticks
forest$get_sticks()
#> [[1]]
#>  [1]   30   46  201  253  544  659  731  758  884  919 1111 1383 1428
#> 
#> [[2]]
#> [1]  1  2  6 10 26 30
#> 

# search for the forest sticks whose corresponding cells have
# birth times 40 time units at most
forest$get_sticks(40)
#> [[1]]
#> [1]  1  2  6 10 26 30
#> 
```
