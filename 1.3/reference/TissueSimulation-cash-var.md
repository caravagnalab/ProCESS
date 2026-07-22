# Building a simulation status variable

This method builds a simulation status variable.

## Arguments

- variable_description:

  The description of the variable to be built. When
  `variable_description` is the string `"Time"`, the elapsed simulation
  time variable is returned. If `variable_description` is set to a
  species name, then the variable representing the cardinality of the
  species is built. Finally, when the parameter is a species name
  followed by `.` and one among `duplications`, `deaths`, or `switches`,
  the variable representing the number of event of the specified type
  occurred since the computation beginning in the species.

## Value

A variable representing the simulation quantity according to the
parameter `variable_description`.

## Details

This method builds a logic variable representing one of the simulation
quantities among:

- cardinality of a species

- number of event among duplications, deaths, and epigenetic switches

- elapsed evolution time

## See also

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

# get the variable representing the simulation time
sim$var("Time")
#> Time

# get the variable representing the cardinality of A[E1]
sim$var("A[E1]")
#> |A[E1]|

# get the variable representing the cardinality of A[E2]
sim$var("A[E2]")
#> |A[E2]|

# get the variable representing the number of duplications
# from A[E2]
sim$var("A[E2].duplications")
#> |A[E2].duplications|

# get the variable representing the number of epigenetic
# switches from A[E1]
sim$var("A[E1].switches")
#> |A[E1].switches|
```
