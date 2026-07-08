# Setting the germline subject

This method sets the germline subject.

## Value

Set the germline subject.

## Details

The subject must be one among those reported by
[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md).

## See also

[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md)
to get the available germline subjects;
[`MutationEngine$get_active_germline()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_active_germline.md)
to get the active germline subject.

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                                                        

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                                                                       

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                                                             

 [█████████████---------------------------] 30% [00m:01s] Loading RS index                                                                                                            

 [█████████████████████████---------------] 61% [00m:02s] Loading RS index                                                                                                            

 [█████████████████████████████████████---] 90% [00m:03s] Loading RS index                                                                                                            

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            


# set the active germline subject data frame
m_engine$set_germline_subject("NA18941")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            
```
