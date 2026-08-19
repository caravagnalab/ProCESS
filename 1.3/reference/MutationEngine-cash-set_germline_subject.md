# Setting the germline subject

This method sets the germline subject.

## Arguments

- subject:

  The germline subject. It must be one among those reported by
  [`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md).

- quiet:

  An optional Boolean flag to avoid the progress bar (default: `FALSE`).

## Details

The subject must be one among those reported by
[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md).

## See also

[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md),
[`MutationEngine$get_active_germline()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_active_germline.md),
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# set the active germline subject data frame
m_engine$set_germline_subject("NA18941")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    
```
