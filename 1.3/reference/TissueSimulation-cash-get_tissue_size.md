# Getting the simulated tissue size

This method returns the size of the simulated tissue.

## Value

The vector `c(x_size, y_size)` of the simulated tissue.

## See also

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation having size 1200x900
sim <- TissueSimulation(width = 1200, height = 900)

# get the tissue size, i.e., expecting c(1200,900)
sim$get_tissue_size()
#> [1] 1200  900
```
