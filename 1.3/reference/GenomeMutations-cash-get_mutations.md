# Getting genome mutations

This method returns a data frame representing the SIDs in the genome.

## Arguments

- only_leaves:

  A Boolean flag to iterate only over forest leaves (default: `FALSE`).

## Value

A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`,
`alt`, `causes`, and `classes`. Each row represent a SID.
