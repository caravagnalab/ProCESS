# Getting the CNA source allele

This method returns the identifier of the allele from which the CNA is
copied.

## Value

The allele from which CNA is copied.

## Examples

``` r
# create an amplification CNA
amp_cna <- Amplification("X", 20002, 100, 1, 0)

# get allele from which `amp_cna` is copied (i.e., 0)
amp_cna$get_src_allele()
#> [1] 0

# create a deletion CNA
del_cna <- Deletion("Y", 40020, 200, 0)

# the deletions have no sources and the method returns NA
del_cna$get_src_allele()
#> [1] NA
```
