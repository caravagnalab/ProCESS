# Getting forest nodes

This method builds a data frame containing forest nodes.

## Value

A data frame reporting, for each node in the forest, the identified
(column `cell_id`), the ancestor identifier (column `ancestor`), the
node's depth (column `depth`), the name of the sample containing the
node, (column `sample`), and the mutant (column `mutant`), the birth
time (column `birth_time`). Whenever, the simulation has epigenetic
states, the data frame also contains the column `epistate`.

## See also

[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a forest example
forest <- example("SampleForest")

# get the data frame of the nodes
nodes <- forest$get_nodes()

# print the first lines of the data frame
head(nodes)
#>   cell_id ancestor depth mutant epistate sample birth_time
#> 1       1       NA     0      A       E1   <NA>   0.000000
#> 2       2        1     1      A       E1   <NA>   5.741436
#> 3       3        1     1      A       E1   <NA>   5.741436
#> 4       4        3     2      A       E1   <NA>   9.232361
#> 5       5        3     2      A       E1   <NA>   9.232361
#> 6       6        2     2      A       E2   <NA>  15.850958
```
