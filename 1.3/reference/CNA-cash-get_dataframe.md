# Getting the CNA data frame

This method builds a data frame representing the CNA.

## Details

The data frame contains the columns `chr`, `from`, `length`, `alt_base`,
`allele`", `src_allele`, and `type`.

## Examples

``` r
# create an amplification CNA
amp_cna <- Amplification("X", 20002, 100)

amp_cna$get_dataframe()
#>   chr  from length allele src_allele type
#> 1   X 20002    100     NA         NA    A

# create a deletion CNA
del_cna <- Deletion("Y", 40020, 200, 0)

del_cna$get_dataframe()
#>   chr  from length allele src_allele type
#> 1   Y 40020    200      0         NA    D
```
