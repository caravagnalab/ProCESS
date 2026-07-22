# The fragment sequence

This property is the fragment sequence.

## Value

The fragment sequence.

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

# get the fragment sequence
fragment$sequence
#> [1] "TGGCTCACTGCAAGCTCCGCGTCCCGGGTTCATGCCATTCTCCTGCCTCAGGTAGCTGGGACTACTGGCACCCACCACTACGCCCGGATACTTTTTGTAT"
```
