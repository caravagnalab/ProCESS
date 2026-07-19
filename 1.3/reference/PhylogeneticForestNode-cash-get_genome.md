# Getting the corresponding cell genome

This method computes the corresponding cell genome.

## Value

The genome of the corresponding cell.

## Details

This method computes the corresponding cell genome by browsing the whole
forest branch from the root down to the node. Whenever the genomes of
many nodes are needed, using the node tour with genomes is preferable.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# the genome only has chromosome 22 because the forest was
# built by using the setup "demo"
node$get_genome()
#> GenomeMutations: 1 chrs 6 alleles
```
