# Simulating cell evolution

This method simulates cell evolution until a formula does not hold.

## Arguments

- formula:

  The formula that will be satisfied at the end of the simulation.

- quiet:

  An optional Boolean flag to avoid the progress bar (default: `FALSE`).

## See also

[`TissueSimulation$var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md),
[`TissueSimulation$run_up_to_time()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_time.md),
[`TissueSimulation$run_up_to_event()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_event.md),
[`TissueSimulation$run_up_to_size()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_size.md),
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

# get the variable representing the simulation time
v_time <- sim$var("Time")

# get the variable representing the cardinality of A[E1]
va_e1 <- sim$var("A[E1]")

# get the variable representing the cardinality of A[E2]
va_e2 <- sim$var("A[E2]")

# get the variable representing the number of epigenetic
# switches from A[E1]
va_ps <- sim$var("A[E1].switches")

# build a condition stating that the cardinality of A[E1]
# doubles that of A[E2]
c1 <- va_e1 >= 2*va_e2

# build a condition that holds when there are more than
# 100000 live cells of mutant A
c2 <- va_e1 + va_e2 > 1e5

# build a condition that holds when less than 4000 switched
# from A[E1] have occurred
c3 <- va_ps < 4000

# build a condition that holds when 40 time unit have been
# simulated at least
c4 <- v_time >= 40

# build a condition that holds when c4 and at least one
# among c1, c2, and c3 hold
c5 <- c4 & (c1 | c2 | c3)
c5
#> Time>=40 and (|A[E1]|>=2*|A[E2]| or |A[E1]|+|A[E2]|>100000 or 4000>|A[E1].switches|)

# run the simulation while c5 does not hold
sim$run_until(c5)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


sim
#> ──  ProCESS   D   S   M  ProCESS_20260811-224500 ───── ▣  [1000x1000]  ⏱ 40.1 ──
#> 
#> ── Species: 2, with epigenetics 
#>    
#>    =======  ====  ====  ======  =========
#>    species   λ     δ    counts      %    
#>    =======  ====  ====  ======  =========
#>      A[E1]  0.20  0.01   528    91.986063
#>      A[E2]  0.08  0.01    46    8.013937 
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
#> ── Firings: 723 total 
#> 
#>  Species A[E1]:  48 (deaths), 609 (duplications) and  36 (switches)
#>  Species A[E2]:   8 (deaths),  20 (duplications) and   2 (switches)
#> ✖ The simulation has no samples yet!
```
