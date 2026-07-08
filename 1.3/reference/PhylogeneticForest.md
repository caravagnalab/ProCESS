# The phylogenetic forest of cells in samples

This class represents the phylogenetic forest of the cells sampled
during the computation.

## Details

The leaves of his forest are the sampled cells. This class is analogous
to the class
[`SampleForest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
but each node is labelled with the mutations occurring for the first
time on the cell represented by the node itself. Moreover each leaf is
also associated with the genome mutations occurring in the corresponding
cell.
