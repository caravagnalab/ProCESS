# Getting the simulation name

This method returns the simulation name

## Value

The simulation name, which corresponds to the name of the directory in
which the simulation is saving its progresses.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# Expecting "test"
sim$get_name()
#> [1] "ProCESS_20260811-224433"
```
