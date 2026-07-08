# Represent a simulation quantity

The objects of this class represent one among the following quantities:

- the cardinality of a species;

- the number fired event among deaths, duplications and switches in a
  species;

- the elapse simulation time.

## See also

[`TissueSimulation$var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md)
to build a variable

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

# get the variable representing the cardinality of "A[E1]" in sim
sim$var("A[E1]")
#> |A[E1]|

# get the variable representing the number of duplications in "A[E1]"
sim$var("A[E1].duplications")
#> |A[E1].duplications|

# get the variable representing the simulation time
sim$var("Time")
#> Time

# the logic variables can be stored in an R variable
va_p <- sim$var("A[E1]")
va_p
#> |A[E1]|
```
