# Simulating cell evolution

This method simulates cell evolution until the number of cells in a
species reaches a specified threshold.

## Arguments

- species:

  The species whose number of cells is considered.

- num_of_cells:

  The threshold for the cell number.

- quiet:

  An optional Boolean flag to avoid the progress bar (default: FALSE).

## See also

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md),
[`TissueSimulation$run_up_to_time()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_time.md),
[`TissueSimulation$run_up_to_event()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_event.md),
[`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)

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

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# simulate the tissue until the species "A[E2]" account for 100
# contemporary cells
sim$run_up_to_size(species = "A[E2]", num_of_cells = 100)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


sim$get_counts()
#>   mutant epistate counts overall
#> 1      A       E1   1181    2843
#> 2      A       E2    100     176
```
