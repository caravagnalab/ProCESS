# Getting the bulk allelic fragmentation data frame

This method returns a data frame representing the bulk allelic
fragmentation of the genome.

## Arguments

- sample_name:

  The name of the sample whose bulk allelic fragmentation is aimed
  (optional).

## Value

A data frame reporting, for each genomic fragment and for all the
allelic type on the genomic fragment, the chromosome (`chr`), the first
base position (`begin`), the last base position (`end`), the number of
copy of the major and minor alleles (`major` and `minor`, respectively),
and the ratio between the number of cells exhibiting this allelic type
and the total number of cells in the sample.

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)
for usage examples.
