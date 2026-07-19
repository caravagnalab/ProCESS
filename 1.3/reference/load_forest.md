# Loading forests

This method loads a forest from a file.

## Usage

``` r
load_forest(filename, quiet)
```

## Arguments

- filename:

  The path of the file from which the forest must be loaded.

- quiet:

  An optional Boolean flag to avoid the progress bar (default: FALSE).

## Value

The loaded forest

## Details

This method loads a forest, being either a sample forest or a
phylogenetic forest, from a file.

## See also

[`SampleForest$save()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-save.md),
[`PhylogeneticForest$save()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-save.md),
[`load_sample_forest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_sample_forest.md),
[`load_phylogenetic_forest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_phylogenetic_forest.md)
