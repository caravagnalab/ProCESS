# Computing the forest sticks

This method computes the forest sticks.

## Arguments

- birth_threshold:

  The maximum birth time for the cells associated to the returned sticks
  (optional).

## Value

The list of the forest sticks whose associated cells have birth time
smaller than or equal to `birth_threshold`. Each stick is represented as
the list of cell identifiers labelling the nodes in the stick from the
higher to the deeper in the forest.

## Details

A *crucial node* of a forest is a root of the forest, a node whose
parent belongs to a different species, or the most recent common
ancestor of two crucial nodes.

A *stick* is a path of the forest in which the only crucial nodes are
the first and the last one.

This method returns the list of the forest sticks. Each stick is
represented by the sequence of cell identifiers labelling the nodes in
the stick.

## See also

[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md)

## Examples

``` r
# use a forest example
forest <- example("SampleForest")

# search for the forest sticks
head(forest$get_sticks())
#> [[1]]
#> [1] 156 306
#> 
#> [[2]]
#> [1] 102 156
#> 
#> [[3]]
#> [1] 1848 1924
#> 
#> [[4]]
#> [1] 18904 19052
#> 
#> [[5]]
#>  [1]   6615   7064   7274   7282   7606   7855   8339   8406  10019  12170
#> [11]  15672  55124  68994 103049 109325 118301 187675 196109 206870
#> 
#> [[6]]
#> [1]  6615  7065  7184  9859 11810 16729 17848 18743 18904
#> 

# search for the forest sticks whose first node corresponding cells have
# birth times 40 time units at most
forest$get_sticks(40)
#> [[1]]
#> [1]  8 10 13 21 36 49 68
#> 
#> [[2]]
#> [1]  8 11 17 22
#> 
#> [[3]]
#> [1] 6 8
#> 
#> [[4]]
#> [1] 1 2 6
#> 
#> [[5]]
#> [1]  1  3  5 18 34
#> 
```
