# Getting genome mutations

This method returns a data frame representing the SIDs in the genome.

## Arguments

- with_germline:

  A Boolean flag to add germline mutations (default: `TRUE`).

## Value

A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`,
`alt`, `cause`, and `nature`. Each row represent a SID.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get the first mutations in the genome
head(genome$get_mutations())
#>   chr     from allele            ref  alt cause         nature
#> 1  22 16085675      0      GCCTCCCGA    G     A         driver
#> 2  22 16095655      0              T    A  SBS1 pre-neoplastic
#> 3  22 16099091      0              G GATA   ID1 pre-neoplastic
#> 4  22 16165397      0              T  TTG   ID1 pre-neoplastic
#> 5  22 16241450      0              T    G  SBS1 pre-neoplastic
#> 6  22 16288596      0 CGCGTGCGGCGTGC    C   ID1 pre-neoplastic
```
