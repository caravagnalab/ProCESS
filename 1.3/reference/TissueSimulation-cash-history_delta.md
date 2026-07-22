# The delta time between time series samples

This value is the maximum time between two successive time series data
samples.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# get the delta time between two time series samples (0 by default)
sim$history_delta
#> [1] 0

# set the delta time between two time series samples
sim$history_delta <- 20
```
