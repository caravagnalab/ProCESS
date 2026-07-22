# Getting the sampled cells' CNAs

This method returns the CNAs of the sample cells.

## Value

A data frame reporting `cell_id`, `type` (`"A"` for amplifications and
`"D"` for deletions), `chr`, `begin` (i.e., the first CNA locus in the
chromosome), `end` (i.e., last CNA locus in the chromosome), `allele`,
`src allele` (the allele origin for amplifications, `NA` for deletions),
and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
`"pre-neoplastic"`).

## Details

This method builds a data frame representing all the CNAs in the cells
sampled during the simulation and represented by the leaves of the
phylogenetic forest.

## See also

[`PhylogeneticForest$get_sampled_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the sampled cell CNAs
CNAs <- forest$get_sampled_cell_CNAs()

# print the first lines of the data frame
head(CNAs)
#>   chr    begin      end type allele src.allele cause nature cell_id
#> 1  22 10303470 10503469    A      2          0       driver  183352
#> 2  22  5010000  5209999    D      1         NA       driver  183352
#> 3  22 10303470 10503469    A      2          0       driver  165816
#> 4  22  5010000  5209999    D      1         NA       driver  165816
#> 5  22 10303470 10503469    A      2          0       driver  177858
#> 6  22  5010000  5209999    D      1         NA       driver  177858
```
