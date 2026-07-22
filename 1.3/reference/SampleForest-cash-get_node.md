# Getting a node of the forest

This method returns the node of the forest

## Arguments

- cell_id:

  The identifier of the cell whose node is aimed.

## Value

The
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
object associated to the cell whose identifier is `cell_id`.

## Details

This method returns the node of the forest whose corresponding cell has
a specified identifier.

## See also

[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a forest example
forest <- example("SampleForest")

# get the node corresponding to the cell having 2 as cell identifier
forest$get_node(2)
#> SampleForestNode(cell_id = 2, species = "A[E1]")
```
