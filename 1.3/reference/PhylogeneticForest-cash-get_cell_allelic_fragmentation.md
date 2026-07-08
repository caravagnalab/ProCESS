# Getting the cell allelic fragmentation data frame

This method returns a data frame representing the allelic fragmentation
of each sampled cell.

## Value

A data frame reporting, for each cell, for each genomic fragment, and
for all the allelic type on the genomic fragment, the cell identifier
(`cell_id`), the chromosome (`chr`), the first base position (`begin`),
the last base position (`end`), and the number of copy of the major and
minor alleles (`major` and `minor`, respectively).

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)
for usage examples.
