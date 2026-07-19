# Getting genome mutations

This method returns a data frame representing the CNAs in the genome.

## Value

A data frame consisting of 8 columns: `chr`, `begin`, `end`, `type`,
`allele`, `src.allele`, `cause`, and `nature`. Each row represent a SID.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get the first CNAs in the genome
head(genome$get_CNAs())
#>   chr    begin      end type allele src.allele cause nature
#> 1  22 10303470 10503469    A      2          0       driver
#> 2  22  5010000  5209999    D      1         NA       driver
```
