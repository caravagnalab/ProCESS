# Getting the mutations laying on a fragment

This method returns the data frame of the mutations laying on the
fragment.

## Value

A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`,
`alt`, `cause`, and `nature`. Each row represent a SID.

## See also

[`GenomeFragment`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get a fragment
fragment <- genome$get_fragment("22", 0, 16085625, 100)

# get the fragment mutations
fragment$get_mutations()
#>   chr allele     from       ref alt cause nature
#> 1  22      0 16085675 GCCTCCCGA   G     A driver
```
