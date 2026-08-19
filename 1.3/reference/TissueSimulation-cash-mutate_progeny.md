# Generating a mutated progeny

This method generates a mutated progeny.

## Arguments

- cell_position:

  The position of the cell whose offspring will mutate.

- mutated_mutant:

  The mutant of the mutated cell.

## Details

It simulates both the duplication of the cell in the specified position
and the birth of one cells of a given mutant that preserves the
epigenetic status of the original cell. The mutated cell will be located
in the position of its parent.

## See also

[`TissueSimulation$choose_cell_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-choose_cell_in.md),
[`TissueSimulation$choose_border_cell_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-choose_border_cell_in.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# add mutant "B" and set its species rates
sim$add_mutant("B",
               list(E1 = list(duplication = 0.15, death = 0.3, E2 = 0.1),
                    E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                    


# get the number of cells per species. No cell in "B" yet.
sim$get_counts()
#>   mutant epistate counts overall
#> 1      A       E1   2043    5150
#> 2      A       E2    184     340
#> 3      B       E1      0       0
#> 4      B       E2      0       0

# duplicate the cell in position (503, 492). One of
# its direct descendants will have mutant "B"
# sim$mutate_progeny(503, 492, "B")

# the output of `choose_cell_in`, `choose_border_cell_in` and `get_cell`
# can also be used as input for `mutate_progeny`
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

# get the number of cells per species again.
# Now, "B" consists of one cell
sim$get_counts()
#>   mutant epistate counts overall
#> 1      A       E1   2043    5152
#> 2      A       E2    184     340
#> 3      B       E1      1       1
#> 4      B       E2      0       0
```
