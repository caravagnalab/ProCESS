# Getting the sampled cell mutations

The mutations of the sampled cell.

## Details

This property contains a data frame that represents the sampled cell
mutations. The data frame format is analogous to that returned by
`PhylogeneticForest$get_sampled_cell_mutations()`: it has columns
`cell_id`, `chr`, (i.e., the mutation chromosome), `from` (i.e.,
position in the chromosome), `allele` (in which the mutation occurs),
`ref`, `alt`, `type` (i.e., either `"SNV"` or `"indel"`), `cause`, and
`class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
`"pre_neoplastic"`).

## See also

[`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md),
and
[`vignette("sample_partition")`](https://caravagnalab.github.io/ProCESS/1.3/articles/sample_partition.md).
