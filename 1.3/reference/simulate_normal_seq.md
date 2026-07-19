# Simulating wild-type sequencing

This method simulates a wild-type sample sequencing in a phylogenetic
forest. Add the cells in the wild-type sample contains the germline
mutations. The forest pre-neoplastic mutations are also added to the
sample by default. However, they can be avoided by using the parameter
`with_pre_neoplastic`.

## Arguments

- phylo_forest:

  A phylogenetic forest.

- sequencer:

  The sequencer that performs the sequencing simulation (default: an
  `ErrorlessIlluminaSequencer`).

- reference_genome:

  The reference genome (default: NULL to use the mutation engine
  reference genome).

- chromosomes:

  The chromosomes that must be considered (default: `NULL`, i.e., all
  the reference chromosomes).

- coverage:

  The sequencing coverage (default: \\10\\).

- read_size:

  The read size (default: \\150\\).

- insert_size_mean:

  The insert size mean. Use 0 for single read sequencing and any value
  greater than 0 for pair read sequencing (default: \\0\\).

- insert_size_stddev:

  The insert size standard deviation. (default: \\10\\).

- output_dir:

  The SAM output directory (default: `"ProCESS_normal_SAM"`).

- write_SAM:

  A Boolean flag to enable/disable SAM generation (default: `TRUE`).

- update_SAM:

  Update the output directory (default: `FALSE`).

- filename_prefix:

  The prefix of the output SAM file name (default: `"chr_"`).

- template_name_prefix:

  The template name prefix (default: `"r"`).

- missed_SID_statistics:

  A Boolean flag to collect statistics also about the mutations that are
  not covered by any of the simulated reads (default: `FALSE`).

- germline_statistics:

  A Boolean flag to collect statistics also about the germinal mutations
  that are not covered by any of the simulated reads (default: `FALSE`).

- wide_format:

  A Boolean flag to request wide/long format for the mutation output
  (default: `TRUE`).

- seed:

  The random seed for the internal random generator (optional).

- quiet:

  A Boolean flag to enable/disable the progress bar (default: FALSE).

## Value

A named list of two elements: the sequencing output data frame (name
`mutations`) and the calling parameters (name `parameters`).

If `wide_format` is set to `true`, the sequencing output data frame
reports, for each of the observed SNVs and indels, the chromosome and
the position in which it occurs (columns `chr` and `from`), the
reference and the alternative sequence, the cause, and the nature of the
mutation (columns `ref`, `alt`, `cause`, and `nature`, respectively).
Moreover, the returned data frame contains three columns: the number of
reads in which the corresponding SNV occurs (column `normal.sample.NV`),
the coverage of the SNV locus (column `normal.sample.DP`), and the
corresponding VAF (column `normal.sample.VAF`).

Instead, when `wide_format` is set to `false`, the output data frame
contains a row for each mutation in each sample and consists of 10
columns: `sample`, `chr`, `from`, `ref`, `alt`, `cause`, `nature`, `NV`,
`DP`, and `VAF`. The column `sample` contains the name of the sample in
which the mutation has been identified. The columns `chr`,
from`, `ref`, `alt`, `cause`, and `nature` correspond to those of the wide_format`
output. The columns `NV`, `DP`, and `VAF` maintain the number of
occurrences, the coverage, and the VAF of the mutation in cited sample.

## See also

[`BasicIlluminaSequencer()`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer.md)
and
[`ErrorlessIlluminaSequencer()`](https://caravagnalab.github.io/ProCESS/1.3/reference/ErrorlessIlluminaSequencer.md)
as sequencer types, and
[`vignette("sequencing")`](https://caravagnalab.github.io/ProCESS/1.3/articles/sequencing.md)
for usage examples.
