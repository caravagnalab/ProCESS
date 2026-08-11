# Switching on and off the infinite sites model

This property enables/disables the infinite sites model.

## Details

When it is `TRUE`, the infinite sites model is enabled and new mutations
are exclusively placed in mutation-free loci.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# create a demonstrative mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# the infinite sites model is enabled by default
m_engine$infinite_sites_model
#> [1] TRUE

# the infinite sites model can be disabled
m_engine$infinite_sites_model <- FALSE

m_engine$infinite_sites_model
#> [1] FALSE
```
