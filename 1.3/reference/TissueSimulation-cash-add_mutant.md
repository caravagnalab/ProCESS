# Adding a mutant and its species

This method adds a mutant and its species to the simulation.

## Arguments

- mutant:

  The mutant name.

- rate_list:

  The list of the mutant's rates (optional).

## Details

This method adds a mutant to the simulation. The method also creates the
species associated to the new mutant according to the known epigenetic
states. The default rate of the new species is set to zero. Optionally,
user can provide a list specifying the rates of the associated species.

## See also

[`TissueSimulation$add_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutants.md),
[`TissueSimulation$add_epigenetic_state()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_state.md),
[`TissueSimulation$add_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_states.md),
[`TissueSimulation$get_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_mutants.md),
[`TissueSimulation$set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md),
[`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# see the simulation setup
sim
#> ──  ProCESS   D   S   M  ProCESS_20260708-171444 ────────────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ✖ The simulation has no samples yet!

# add the mutant "A" to the simulation.
sim$add_mutant(name = "A")

# see the simulation setup
sim
#> ──  ProCESS   D   S   M  ProCESS_20260708-171444 ────────────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ── Species: 1, without epigenetics 
#>    
#>    =======  ===  ===  ======  ===
#>    species   λ    δ   counts   % 
#>    =======  ===  ===  ======  ===
#>          A   0    0     0     NaN
#>    =======  ===  ===  ======  ===
#> 
#> ── Firings: 0 total 
#> 
#> ✖ The simulation has no samples yet!

# add the mutant "B" to the simulation and set its duplication and death rates
sim$add_mutant(name = "B", rate_list = c(duplication = 0.3, death = 0.1))

# see the simulation setup
sim
#> ──  ProCESS   D   S   M  ProCESS_20260708-171444 ────────────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ── Species: 2, without epigenetics 
#>    
#>    =======  ===  ===  ======  ===
#>    species   λ    δ   counts   % 
#>    =======  ===  ===  ======  ===
#>          A  0.0  0.0    0     NaN
#>          B  0.3  0.1    0     NaN
#>    =======  ===  ===  ======  ===
#> 
#> ── Firings: 0 total 
#> 
#> ✖ The simulation has no samples yet!

# add epigenetic states (rates are reset)
sim$add_epigenetic_states(c("E1", "E2", "E3"))

# add the mutant "C" to the simulation, set the duplication and death rates
# of all its species, and differentiate "C[E1]" by setting its death rate
# and the rates of the switch toward "C[E2]" and "C[E3]".
sim$add_mutant("C", list(duplication = 0.3, death = 0.1,
                         E1=list(death = 0.2, E2=0.01, E3=0.1)))

# see the simulation setup
sim
#> ──  ProCESS   D   S   M  ProCESS_20260708-171444 ────────────────────────────────────────────────────────────────────────────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ── Species: 9, with epigenetics 
#>    
#>    =======  ===  ===  ======  ===
#>    species   λ    δ   counts   % 
#>    =======  ===  ===  ======  ===
#>      A[E1]  0.0  0.0    0     NaN
#>      A[E2]  0.0  0.0    0     NaN
#>      A[E3]  0.0  0.0    0     NaN
#>      B[E1]  0.0  0.0    0     NaN
#>      B[E2]  0.0  0.0    0     NaN
#>      B[E3]  0.0  0.0    0     NaN
#>      C[E1]  0.3  0.2    0     NaN
#>      C[E2]  0.3  0.1    0     NaN
#>      C[E3]  0.3  0.1    0     NaN
#>    =======  ===  ===  ======  ===
#> 
#> ── Epigenetic switches 
#>    
#>    =======  ====  =====
#>    species   ε     dest
#>    =======  ====  =====
#>      C[E1]  0.01  C[E2]
#>      C[E1]  0.10  C[E3]
#>    =======  ====  =====
#> 
#> ── Firings: 0 total 
#> 
#> ✖ The simulation has no samples yet!
```
