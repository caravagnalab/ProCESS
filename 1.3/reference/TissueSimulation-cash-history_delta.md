# The delta time between time series samples

The maximum time between two consecutive time series data samples.

## Details

This property is the maximum time between two consecutive time series
data samples. Notice that this property differs from
[`TissueSimulation$snapshot_triggers`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_triggers.md).

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md),
[`TissueSimulation$snapshot_triggers`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_triggers.md)

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
