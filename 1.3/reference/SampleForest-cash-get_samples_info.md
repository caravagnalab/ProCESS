# Retrieving the samples' information

This method retrieves information about the samples whose cells were
used as leaves of the sample forest.

## Value

A data frame containing, for each sample collected during the
simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`, `ymax`,
`xmax`, `tumour_cells`, and `tumour_cells_in_bbox`. The columns `ymin`,
`xmin`, `ymax`, `xmax` report the boundaries of the sample bounding box,
while `tumour_cells` and `tumour_cells_in_bbox` are the number of tumour
cells in the sample and in the bounding box, respectively.

## See also

[`PhylogeneticForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_samples_info.md)
for usage examples,
[`TissueSimulation$sample_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-sample_cells.md),
[`TissueSimulation$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get information about the samples whose cells
# are the forest leaves
forest$get_samples_info()
#>    name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 S_1_1  0  480  480  490  490          120                  120 552.3795
#> 2 S_1_2  1  500  500  510  510          119                  119 552.3795
#> 3 S_2_1  2  370  551  394  575          613                  613 741.0625
#> 4 S_2_2  3  420  276  444  300          572                  572 741.0625
```
