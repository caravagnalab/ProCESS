# Getting the name of the corresponding cell mutant

This property is the name of the corresponding cell mutant.

## Value

The name of the corresponding cell epigenetic state.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the name of the corresponding cell species
node$species_name
#> [1] "A[E1]"

# get the corresponding cell mutant
node$mutant_name
#> [1] "A"
```
