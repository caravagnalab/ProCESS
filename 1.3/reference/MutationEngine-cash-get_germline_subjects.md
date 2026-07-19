# Getting the germline subjects

This method returns the available germline subjects.

## Value

A data frame the available germline subjects.

## Details

The germline subjects method returns a data frame containing the
available germline subjects. The column `sample` reports the subject
name; the columns `pop` and `super_pop` contain the subject population
and super population, respectively; the column `gender` declares the
subject gender.

## See also

[`MutationEngine$get_active_germline()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_active_germline.md)
to get the available germline subjects;
[`MutationEngine$set_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-set_germline_subject.md)
to set the active germline.

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index         

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded        

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index              

 [█████████████---------------------------] 32% [00m:01s] Loading RS index             

 [██████████████████████████--------------] 64% [00m:02s] Loading RS index             

 [██████████████████████████████████████--] 93% [00m:03s] Loading RS index             

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded             

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline              

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded             


# get the active germline subject data frame
head(m_engine$get_germline_subjects(), 5)
#>    sample pop super_pop gender
#> 1 NA18941 JPT       EAS female
#> 2 NA20513 TSI       EUR   male
```
