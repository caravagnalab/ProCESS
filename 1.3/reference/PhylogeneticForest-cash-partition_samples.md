# Partitioning forest samples

This method partitions the samples in a phylogenetic forest.

## Arguments

- labelling_function:

  An R labelling function that maps any sampled cell to a labelling
  string.

## Details

This method partitions the samples in a phylogenetic forest according to
a labelling function. It works in-place by altering the phylogenetic
forest from which the method is called.

## See also

[`SampleForest$get_subforest_for()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_subforest_for.md)
for usage examples.
