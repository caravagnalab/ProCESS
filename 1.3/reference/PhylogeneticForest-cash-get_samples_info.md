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

[`SampleForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_samples_info.md)
for usage examples,
[`TissueSimulation$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)
