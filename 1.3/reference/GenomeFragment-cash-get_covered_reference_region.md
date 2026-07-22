# Getting the reference region covered by the fragment

This method returns the reference region covered by the fragment.

## Value

A named list representing the reference region covered by the fragment.
The names of the list are `chr`, `from`, and `size`. They stores the
name of the chromosome from which the fragment comes, the position in
the chromosome of the fragment first base, and the size of the covered
region, respectively.

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

# get the covered reference region
fragment$get_covered_reference_region()
#> $chr
#> [1] "22"
#> 
#> $allele
#> [1] 0
#> 
#> $from
#> [1] 16085625
#> 
#> $size
#> [1] 108
#> 
```
