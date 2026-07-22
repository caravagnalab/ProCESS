# Getting the cell allelic fragmentation data frame

This method returns a data frame representing the allelic fragmentation
of each sampled cell.

## Value

A data frame reporting, for each cell, for each genomic fragment, and
for all the allelic type on the genomic fragment, the cell identifier
(`cell_id`), the chromosome (`chr`), the first base position (`begin`),
the last base position (`end`), and the number of copy of the major and
minor alleles (`major` and `minor`, respectively).

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# print the first rows of the cell allelic fragmentation
head(forest$get_cell_allelic_fragmentation())
#>   cell_id chr    begin      end major minor
#> 1  183352  22        1  5009999     2     2
#> 2  183352  22  5010000  5209999     2     0
#> 3  183352  22  5210000 10303469     2     2
#> 4  183352  22 10303470 10503469     4     2
#> 5  183352  22 10503470 51304566     2     2
#> 6  165816  22        1  5009999     2     2
```
