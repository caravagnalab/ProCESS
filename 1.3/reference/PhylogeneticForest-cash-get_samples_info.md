# Retrieving sample information

This method retrieves information about the samples whose cells were
used as leaves of the sample forest.

## Value

A data frame containing, for each sample collected during the
simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`, `ymax`,
`ymax`, `xmax`, `tumour_cells`, `tumour_cells_in_bbox`, `DNA_quantity`,
and `equivalent_normal_cells`. The columns `ymin`, `xmin`, `ymax`, and
`xmax` report the boundaries of the sample bounding box, while
`tumour_cells` and `tumour_cells_in_bbox` are the number of tumour cells
in the sample and in the bounding box, respectively. Finally,
`DNA_quantity` contains the overall number of tumoral bases, i.e., the
sum of the lengths of all the alleles of all the sample tumoral cells,
and `equivalent_normal_cells` contains the number of normal cells that
contain as much DNA as the sample tumoral cells. The `DNA_quantity` is
stored as a real number despite being a natural number as it usually
exceeds the largest natural number that can be natively represented by
R.

## See also

[`SampleForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_samples_info.md),
[`TissueSimulation$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)

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
```
