# Getting the active germline subject

This method returns the active germline subject.

## Value

A data frame the active germline subject.

## Details

The active germline subject is returned as a data frame in which the
column `sample` reports the subject name, the columns `pop` and
`super_pop` contain the subject population and super population,
respectively, and the column `gender` declares the subject gender.

## See also

[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md)
to get the available germline subjects;
[`MutationEngine$set_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-set_germline_subject.md)
to set the active germline subject.

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

 [████████████████████████----------------] 58% [00m:02s] Loading RS index                                                                                                            

 [███████████████████████████████████-----] 86% [00m:03s] Loading RS index                                                                                                            

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            


# get the active germline subject data frame
head(m_engine$get_active_germline(), 5)
#>    sample pop super_pop gender
#> 1 NA18941 JPT       EAS female
```
