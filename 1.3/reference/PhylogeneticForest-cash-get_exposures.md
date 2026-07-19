# Getting the timed exposure data frame

This method returns a data frame representing the exposure evolution
over time.

## Value

A data frame reporting `time`, `signature`, `exposure` and, `type`.

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the exposures used to build the forest
forest$get_exposures()
#>    time signature exposure  type
#> 1     0      SBS1      0.8   SNV
#> 2     0     SBS13      0.2   SNV
#> 3   100      SBS2      0.2   SNV
#> 4   100      SBS3      0.5   SNV
#> 5   100      SBS5      0.3   SNV
#> 6   120      SBS2      0.2   SNV
#> 7   120      SBS3      0.5   SNV
#> 8   120      SBS5      0.3   SNV
#> 9     0      ID13      0.2 indel
#> 10    0       ID2      0.6 indel
#> 11    0      ID21      0.2 indel
#> 12  120       ID1      0.8 indel
#> 13  120       ID9      0.2 indel
```
