# Placing one cell in the tissue

This method places a cell in the tissue.

## Arguments

- species:

  The name of the new cell species.

- x:

  The position on the x axis of the cell.

- y:

  The position on the y axis of the cell.

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell of species "A" in position (500,500)
sim$place_cell("A", 500, 500)
# create a simulation
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.1)))

# place a cell of species "A[E1]" in position (500,500)
sim$place_cell("A[E1]", 500, 500)
```
