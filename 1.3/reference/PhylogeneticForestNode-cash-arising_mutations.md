# Getting the mutations arising in the corresponding cell

This property is the data frame of the mutations arising the
corresponding cell.

## Value

The data frame of the mutations arising the corresponding cell.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the mutations arising in the corresponding cell
node$arising_mutations
#>   order type CNA_type chr    start      end ref alt allele src.allele    nature
#> 1     1  SID     <NA>  22 22301874 22301874   A   C      1         NA passenger
#>   cause code
#> 1  SBS1 <NA>
```
