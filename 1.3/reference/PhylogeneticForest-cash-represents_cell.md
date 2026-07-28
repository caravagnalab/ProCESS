# Testing whether a cell is represented by a forest

This method checks whether a cell is represented by a forest.

## Arguments

- cell_id:

  The identifier of the cell whose presence in the forest is tested.

## Value

`TRUE` if and only if the the forest contains a node whose corresponding
cell identifier is `cell_id`.

## Details

This method returns `TRUE` if and only if the forest contains a node
whose corresponding cell has the specified identifier.

## See also

[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a forest example
forest <- example("PhylogeneticForest")

# load dplyr to ease the code
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# test whether the forest represents cell 2
forest$represents_cell(2)
#> [1] TRUE

# get the node 2's information
forest$get_nodes() %>% dplyr::filter(cell_id == 2)
#>   cell_id ancestor depth mutant epistate sample birth_time
#> 1       2        1     1      A       E1   <NA>   5.741436

# test whether the forest represents cell 2e10
forest$represents_cell(2e10)
#> [1] FALSE

# get the node 2e10's information
forest$get_nodes() %>% dplyr::filter(cell_id == 2e10)
#> [1] cell_id    ancestor   depth      mutant     epistate   sample     birth_time
#> <0 rows> (or 0-length row.names)
```
