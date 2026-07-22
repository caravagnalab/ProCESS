# Internal cells duplication

This property switches between homogeneous and border driven growth
models.

## Details

This Boolean flag switches between homogeneous and border driven growth
models. When it is set to `TRUE`, the border-growth model is used.
Otherwise, the homogeneous-growth model is applied. It is set to `TRUE`
by default.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# is the simulation using the border driven growth model
# (default: TRUE)
sim$border_growth_model
#> [1] TRUE

# switch to the homogeneous-growth model
sim$border_growth_model <- FALSE

# now it is set to `FALSE`
sim$border_growth_model
#> [1] FALSE

# switch back to the border-growth model
sim$border_growth_model <- FALSE
```
