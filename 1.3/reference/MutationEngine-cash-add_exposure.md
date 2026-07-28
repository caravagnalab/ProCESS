# Adding an exposure to the mutation engine

This method adds an exposure to the mutation engine.

## Arguments

- time:

  The simulated time at which the exposure is adopted (optional).

- exposure:

  An exposure for the specified mutation type, i.e., a discrete
  probability distribution over a set of signature. The indel and SNV
  exposures can be specified in the same list.

## Details

The exposure will be used to establish the probability for a passenger
mutation to occur depending on its context.

Each exposure is associated to a time that is the simulated time in
which the set is adopted. If a time is provided the exposure is used
from the specified time on up to the successive exposure change. When an
exposure is added to the mutation engine without specifying the time,
its time is 0.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# create a demonstrative mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
#> Downloading signature files...
#> Signature file downloaded
#> Downloading driver mutation file...
#> Driver mutation file downloaded
#> Decompressing driver mutation file...done
#> Downloading passenger CNAs file...
#> Passenger CNAs file downloaded
#> Decompressing passenger CNAs file...done
#> Downloading germline...
#> Germline downloaded
#> Decompressing mutations...
#> done
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                   

 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22                                                                  

 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22                                                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                

#> done
#> Building repeated sequence index...
#> 
 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                           

 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] RS index built                                                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                      

 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                      

 [█████████████████████-------------------] 52% [00m:02s] Saving RS index                                                                     
done
#> 
 [████████████████████████████████████████] 100% [00m:02s] RS index saved                                                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Saving germline                                                                      

 [████████████████████████████████████████] 100% [00m:00s] Germline saved                                                                     


# add a default set of coefficients that will be used from simulated
# time 0 up to the successive coefficient change. The indel and SNV
# exposures can be specified in the same list.
m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7, ID2 = 0.2, ID3 = 0.3,
  ID20 = 0.5))

# add a default set of coefficients that will be used from simulated
# time 3.2 up to the end of the simulation.
m_engine$add_exposure(3.2, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))

m_engine
#> MutationEngine
#>  Passenger rates
#> 
#>  Driver mutations
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#>      [0, 3.2[: {"SBS1": 0.7, "SBS13": 0.3}
#>      [3.2, ∞[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}
#> 
#>    indel Timed Exposures
#>      [0, ∞[: {"ID2": 0.2, "ID20": 0.5, "ID3": 0.3}
#> 
```
