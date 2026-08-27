# Simulating cell evolution

This method simulates cell evolution until the number of events that
have occurred to cells of a species reaches a specified threshold.

## Arguments

- event:

  The considered event, i.e., `growth`, `death`, or `switch`.

- species:

  The species whose event number is considered.

- num_of_events:

  The threshold for the event number.

- quiet:

  An optional Boolean flag to avoid the progress bar (default: `FALSE`).

## See also

[`TissueSimulation$run_up_to_time()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_time.md),
[`TissueSimulation$run_up_to_size()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_size.md),
[`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md),
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

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# simulate the cell evolution until the number of epigenetic events from
# the species "A[E2]" is less than 100.
sim$run_up_to_event(event = "switch", species = "A[E2]",
                    num_of_events = 100)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                    


sim
#> ──  ProCESS   D   S   M  ProCESS_20260826-222740 ──────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 149 ──
#> 
#> ── Species: 2, with epigenetics 
#>    
#>    =======  ====  ====  ======  =========
#>    species   λ     δ    counts      %    
#>    =======  ====  ====  ======  =========
#>      A[E1]  0.20  0.01  11704   92.121212
#>      A[E2]  0.08  0.01   1001   7.878788 
#>    =======  ====  ====  ======  =========
#> 
#> ── Epigenetic switches 
#>    
#>    =======  ====  =====
#>    species   ε     dest
#>    =======  ====  =====
#>      A[E1]  0.01  A[E2]
#>      A[E2]  0.01  A[E1]
#>    =======  ====  =====
#> 
#> ── Firings: 25516 total 
#> 
#>  Species A[E1]:  5438 (deaths), 17909 (duplications) and   868 (switches)
#>  Species A[E2]:   484 (deaths),   717 (duplications) and   100 (switches)
#> ✖ The simulation has no samples yet!
```
