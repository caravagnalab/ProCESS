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

[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md),
[`MutationEngine$set_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-set_germline_subject.md),
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# get the active germline subject data frame
head(m_engine$get_active_germline(), 5)
#>    sample pop super_pop gender
#> 1 NA18941 JPT       EAS female
```
