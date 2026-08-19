# Getting the reference genome path

This method returns the reference genome path.

## Value

The reference genome path.

## See also

[`PhylogeneticForest$set_reference_path()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-set_reference_path.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the reference path
forest$get_reference_path()
#> [1] "/private/var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/RtmpS68UCn/temp_libpath342a5a608ba5/ProCESS/extdata/example_ref.fasta"
```
