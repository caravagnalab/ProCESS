# Getting the germinal mutations

This method returns the forest SNVs and indels.

## Value

A data frame reporting `chr`, `from` (i.e., the position in the
chromosome), `allele` (in which the mutation occurs), `ref`, `alt`,
`cause`, `type` (i.e., either `"SNV"` or `"indel"`) and `class` (i.e.,
`"germinal"`).

## Details

Its builds a data frame representing all the germinal SNVs and indels of
the cells represented in the phylogenetic forest. The data frame also
reports the allele in which the mutations occur to support double
occurrences due to CNAs.

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)
for usage examples.
