# First-order simulation state formulas

The objects of this class describes properties of the simulation status.

A formula is:

- a relation among two expressions (operators `<`, `<=`, `==`, `!=`,
  `>=`, `>`);

- the conjunction of two formulas (operator `&`);

- the disjunction of two formulas (operator `|`).

## See also

[`Variable()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Variable.md),
[`TissueSimulation$var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md),
[`vignette("tissue_simulation")`](https://caravagnalab.github.io/ProCESS/1.3/articles/tissue_simulation.md)

## Examples

``` r
# build a simulation and add two species to it
# set the seed of the random number generator
set.seed(0)

# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add the mutant "A" and set its species rates
sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
                                   E2 = 0.01),
                         E2 = list(duplication = 0.08, death = 0.01,
                                   E1 = 0.01)))

# get a formula that holds when the cardinality of the mutant A
# is greater than 1000
f1 <- sim$var("A[E1]") + sim$var("A[E2]") > 1000

# get a formula that holds when the simulated time is 10 at least
f2 <- sim$var("Time") >= 40

# get a formula that holds when the number of duplications doubles
# the switch from A[E1]
f3 <- sim$var("A[E1].duplications") > 2 * sim$var("A[E1].switches")

# combine above formulas by using Boolean operators `&` and `|`
f1 & (f2 | f3)
#> |A[E1]|+|A[E2]|>1000 and (Time>=40 or |A[E1].duplications|>2*|A[E1].switches|)
```
