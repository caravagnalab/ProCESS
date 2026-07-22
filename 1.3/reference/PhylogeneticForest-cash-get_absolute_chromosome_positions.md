# Getting the absolute chromosome positions

This method returns the absolute chromosome positions.

## Value

A data frame reporting the name (column `chr`), the length (column
`length`), the initial absolute position (column `from`), and the final
absolute position (column `to`) of each chromosome.

## Details

Its builds a data frame reporting the name, the length, and the initial
and final absolute positions of each chromosome in the reference genome.

## See also

[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get absolute chromosome positions. Since this forest example was built
# by using one single chromosome, the resulting data frame contains only
# one line
forest$get_absolute_chromosome_positions()
#>   chr   length from       to
#> 1  22 51304566    1 51304566
```
