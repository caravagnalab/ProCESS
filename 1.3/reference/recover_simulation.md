# Loading a simulation

This method loads a simulation from the disk.

## Usage

``` r
recover_simulation(name)
```

## Arguments

- name:

  The name of the simulation to be recovered.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the random seed for repeatability
set.seed(0)

# create a simulation having name "recover_simulation_test" and
# save its snapshots in a local directory
sim <- TissueSimulation("recover_simulation_test",
                        epigenetic_states = c("E1", "E2"),
                        save_snapshots = TRUE)

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
                    E2 = list(duplication = 0.01, death = 0.005, E1 = 0.01)))

# place a cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# simulate up to time 50
sim$run_up_to_time(50)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                       


# show the simulation
sim
#> ──  ProCESS   D   S   M  recover_simulation_test ──────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 50 ──
#> 
#> ── Species: 2, with epigenetics 
#>    
#>    =======  ====  =====  ======  =========
#>    species   λ      δ    counts      %    
#>    =======  ====  =====  ======  =========
#>      A[E1]  0.20  0.050   360    90.225564
#>      A[E2]  0.01  0.005    39    9.774436 
#>    =======  ====  =====  ======  =========
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
#> ── Firings: 892 total 
#> 
#>  Species A[E1]: 224 (deaths), 621 (duplications) and  39 (switches)
#>  Species A[E2]:   3 (deaths),   4 (duplications) and   1 (switches)
#> ✖ The simulation has no samples yet!

# remove the object sim from the environment
rm(list = c("sim"))

# the object pointed by sim does not exist any more
exists("sim")
#> [1] FALSE

# recover the simulation from the directory "recover_simulation_test"
sim <- recover_simulation("recover_simulation_test")

sim
#> ──  ProCESS   D   S   M  recover_simulation_test ──────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 50 ──
#> 
#> ── Species: 2, with epigenetics 
#>    
#>    =======  ====  =====  ======  =========
#>    species   λ      δ    counts      %    
#>    =======  ====  =====  ======  =========
#>      A[E1]  0.20  0.050   360    90.225564
#>      A[E2]  0.01  0.005    39    9.774436 
#>    =======  ====  =====  ======  =========
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
#> ── Firings: 892 total 
#> 
#>  Species A[E1]: 224 (deaths), 621 (duplications) and  39 (switches)
#>  Species A[E2]:   3 (deaths),   4 (duplications) and   1 (switches)
#> ✖ The simulation has no samples yet!

# delete dump directory
unlink("recover_simulation_test", recursive = TRUE)
```
