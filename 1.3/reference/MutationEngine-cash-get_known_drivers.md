# Getting the known driver mutations

This method returns the known driver mutations.

## Value

A data frame containing the known driver.

## Details

The mutation are returned in a data frame reporting the known driver
mutations together with their types, associated tumours, affected genes,
and code name. The first three columns ("`chr`", `from`, and `to`)
report the mutation chromosome, the initial position and the final
position, respectively. The next three columns ("`ref`", `alt`, and
`mutation_type`) describe the reference sequence, the altered sequence,
and the type of the mutation. The last four columns ("`driver_gene`",
`driver_code`, `driver_CDS`, and `tumour_type`) detail the affected
gene, the driver code, which can be used to specify the mutation when
adding a mutant to the mutation engine, the variant code, and the tumour
type associated to the mutation.

## See also

[`MutationEngine$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-add_mutant.md)

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

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

 [██████████████████████████████████████--] 94% [00m:03s] Loading RS index      

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded      

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline       

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded      


# get the known driver data frame
head(m_engine$get_known_drivers(), 5)
#>   chr     from       to ref alt mutation_type driver_gene driver_code
#> 1   1 11168260 11168260   G   A           SNV        MTOR MTOR L2538F
#> 2   1 11168260 11168260   G   A           SNV        MTOR MTOR L2538F
#> 3   1 11168260 11168260   G   A           SNV        MTOR MTOR L2538F
#> 4   1 11168260 11168260   G   A           SNV        MTOR MTOR L2538F
#> 5   1 11168260 11168260   G   A           SNV        MTOR MTOR L2538F
#>   driver_CDS tumour_type
#> 1  c.7612C>T        BLCA
#> 2  c.7612C>T        BRCA
#> 3  c.7612C>T       CCRCC
#> 4  c.7612C>T       CHRCC
#> 5  c.7612C>T      CLLSLL
```
