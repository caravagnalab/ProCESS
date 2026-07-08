# Getting the simulation name

This method returns the simulation name

## Value

The simulation name, which corresponds to the name of the directory in
which the simulation is saving its progresses.

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# Expecting "test"
sim$get_name()
#> [1] "ProCESS_20260708-171501"
```
