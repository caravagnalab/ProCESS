# Partitioning forest samples

This method partitions the samples in a phylogenetic forest.

## Arguments

- labelling_function:

  An R labelling function that maps any forest node to a labelling
  string.

## Details

This method partitions the samples in a phylogenetic forest according to
a labelling function over the corresponding forest nodes. It works
in-place by altering the phylogenetic forest from which the method is
called.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# show information about samples
forest$get_samples_info()
#>    name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 S_1_1  0  480  480  490  490          120                  120 552.3795
#> 2 S_1_2  1  500  500  510  510          119                  119 552.3795
#> 3 S_2_1  2  370  551  394  575          613                  613 741.0625
#> 4 S_2_2  3  420  276  444  300          572                  572 741.0625
#>   DNA_quantity equivalent_normal_cells
#> 1  24672048948                240.4469
#> 2  24383999806                237.6397
#> 3 125944671782               1227.4217
#> 4 117384847008               1144.0000

# the labelling function parameter has type `PhylogeneticForestNode`
birth_time_labelling <- function(node) {
  if (node$birth_time > 421) {
    return("YOUNG")
  }

  if (node$birth_time > 321) {
    return("MIDDLE_AGED")
  }

  return("OLD")
}

# partition the samples according to the labelling function
forest$partition_samples(birth_time_labelling)

# show the new samples
forest$get_samples_info()
#>                name  id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox
#> 1 S_1_1_MIDDLE_AGED 122  480  480  490  490            7                  120
#> 2       S_1_1_YOUNG 125  480  480  490  490          113                  120
#> 3 S_1_2_MIDDLE_AGED 123  500  500  510  510            6                  119
#> 4         S_1_2_OLD 121  500  500  510  510            1                  119
#> 5       S_1_2_YOUNG 124  500  500  510  510          112                  119
#> 6       S_2_1_YOUNG 126  370  551  394  575          613                  613
#> 7       S_2_2_YOUNG 127  420  276  444  300          572                  572
#>       time DNA_quantity equivalent_normal_cells
#> 1 552.3795   1436527848                 14.0000
#> 2 552.3795  23235521100                226.4469
#> 3 552.3795   1231309584                 12.0000
#> 4 552.3795    205218264                  2.0000
#> 5 552.3795  22947471958                223.6397
#> 6 741.0625 125944671782               1227.4217
#> 7 741.0625 117384847008               1144.0000
```
