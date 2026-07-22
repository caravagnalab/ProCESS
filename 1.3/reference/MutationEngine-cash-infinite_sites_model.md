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
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                                                   

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                                                                  

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                                                        

 [████████████----------------------------] 29% [00m:01s] Loading RS index                                                                                                       

 [████████████████████████----------------] 58% [00m:02s] Loading RS index                                                                                                       

 [███████████████████████████████████-----] 87% [00m:03s] Loading RS index                                                                                                       

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                                                       

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                        

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                       


# the infinite sites model is enabled by default
m_engine$infinite_sites_model
#> [1] TRUE

# the infinite sites model can be disabled
m_engine$infinite_sites_model <- FALSE

m_engine$infinite_sites_model
#> [1] FALSE
```
