# Getting the registered species and their rates

This method returns the registered species and their rates.

## Value

A data frame containing the registered species rates.

## Details

The registered species and their rates during the simulation are
returned in a data frame. The column `species` contains the species
names; the columns `time`, `SNV_rate`, `indel_rate`, and `CNA_rate`
store the time from which rates hold, and the corresponding the SNV,
indel, and CNA rates, respectively.

## See also

[`PhylogeneticForest$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_species_info.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                                                        

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                                                                       

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                                                             

 [█████████████---------------------------] 31% [00m:01s] Loading RS index                                                                                                            

 [█████████████████████████---------------] 62% [00m:02s] Loading RS index                                                                                                            

 [█████████████████████████████████████---] 91% [00m:03s] Loading RS index                                                                                                            

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            


# get the active germline subject data frame
head(m_engine$get_species_info(), 5)
#> [1] species    time       SNV_rate   indel_rate CNA_rate  
#> <0 rows> (or 0-length row.names)
```
