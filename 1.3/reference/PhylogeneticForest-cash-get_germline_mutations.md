# Getting the germinal mutations

This method returns the forest SNVs and indels.

## Value

A data frame reporting `chr`, `from` (i.e., the position in the
chromosome), `allele` (in which the mutation occurs), `ref`, `alt`,
`cause`, `type` (i.e., either `"SNV"` or `"indel"`) and `class` (i.e.,
`"germinal"`).

## Details

Its builds a data frame representing all the germinal SNVs and indels of
the cells represented in the phylogenetic forest. The data frame also
reports the allele in which the mutations occur to support double
occurrences due to CNAs.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the first germline mutations
head(forest$get_germline_mutations())
#>   chr     from allele ref   alt cause   nature
#> 1  22 16051493      0   G     A       germinal
#> 2  22 16052167      0   A AAAAC       germinal
#> 3  22 16053659      0   A     C       germinal
#> 4  22 16054740      0   A     G       germinal
#> 5  22 16055942      0   C     T       germinal
#> 6  22 16058070      0   A     G       germinal
```
