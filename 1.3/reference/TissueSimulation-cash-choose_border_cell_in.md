# Picking one border cell

This method chooses one border cell among those belonging to one of the
specified mutants and species.

## Arguments

- names:

  The names of the mutants or species among which the cell must be
  choosen. Can either be a single name or a list of names.

- lower_corner:

  The lower corner of the rectangular selection (optional).

- upper_corner:

  The upper corner of the rectangular selection (optional).

## Value

A list reporting `cell_id`, `mutant`, `position_x`, `position_y`, and,
when the simulation has epigenetic states, `epistate` of the chosen
cell.

## Details

It randomly chooses one of the cells belonging to either a mutant or a
species that has a wild-type cell in its neighborhood. Optionally, the
lower and upper corners of a tissue rectangular selection can be
provided to obtain one cell in the rectangle.

## See also

[`TissueSimulation$choose_cell_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-choose_cell_in.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.1)))

# add mutant "B" and set its species rates
sim$add_mutant("B",
               list(E1 = list(duplication = 0.15, death = 0.1, E2 = 0.1),
                    E2 = list(duplication = 0.4, death = 0.01, E1 = 0.01)))

# schedule a mutation from "A" to "B"
sim$schedule_mutation("A", "B", 20)

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# set the death activation level
sim$death_activation_level <- 100

# run the simulation until "B[E2]" accounts for less than 1000 cells
sim$run_up_to_size("B[E2]", 1000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# Randomly choose one cell in "B" in the tissue
sim$choose_border_cell_in("B")
#>   cell_id mutant epistate position_x position_y
#> 1    5101      B       E2        472        510

# Randomly choose one cell in "B" in a rectangle
sim$choose_border_cell_in("B", c(500, 500), c(520, 520))
#>   cell_id mutant epistate position_x position_y
#> 1    4813      B       E2        504        519

# Randomly choose one cell in "B[E1]"
sim$choose_border_cell_in("B[E1]")
#>   cell_id mutant epistate position_x position_y
#> 1    4692      B       E1        488        523

# Randomly choose one cell in "B[E1]" and any species in "A"
sim$choose_border_cell_in(c("B[E1]", "A"))
#>   cell_id mutant epistate position_x position_y
#> 1    4363      A       E1        512        487
```
