# Getting the forest nodes

This method returns the nodes of the forest.

## Value

A data frame representing, for each node in the forest, the identified
(column `cell_id`), the ancestor identifier (column `ancestor`), the
node's depth (column `depth`), the name of the sample containing the
node (column `sample`), the mutant (column `mutant`), the birth time
(column `birth_time`), and, whenever the simulation has epigenetic
states, the epigenetic state (column `epistate`).

## See also

[`SampleForest$get_nodes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_nodes.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

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

# use a phylogenetic forest example
forest <- example("PhylogeneticForest - no epistates")

# get the data frame of the nodes
nodes <- forest$get_nodes()

# print the first lines of the data frame
head(nodes)
#>   cell_id ancestor depth mutant sample birth_time
#> 1       1       NA     0      A   <NA>   0.000000
#> 2       2        1     1      A   <NA>   7.176795
#> 3       3        1     1      A   <NA>   7.176795
#> 4       4        2     2      A   <NA>  10.681481
#> 5       5        2     2      A   <NA>  10.681481
#> 6       6        3     2      A   <NA>  17.051085
```
