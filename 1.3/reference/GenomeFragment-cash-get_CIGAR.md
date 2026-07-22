# Getting the fragment CIGAR

This method returns the CIGAR code of the fragment with respect to the
reference genome.

## Value

The CIGAR code of the fragment with respect to the reference genome.

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

# get its CIGAR with respect to the reference genome
fragment$get_CIGAR()
#> [1] "51M8D49M"
```
