# Simulating cell evolution

Simulating cell evolution

## Arguments

- time:

  The final simulation time.

- quiet:

  An optional Boolean flag to avoid the progress bar (default: FALSE).

## See also

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md),
[`TissueSimulation$run_up_to_event()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_event.md),
[`TissueSimulation$run_up_to_size()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_size.md),
[`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation without epigenetic states
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", list(duplication = 0.3, death = 0.01))

# place an "A" cell in the tissue
sim$place_cell("A", 500, 500)

# simulate the tissue up to simulate timed 40
sim$run_up_to_time(40)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot             
```
