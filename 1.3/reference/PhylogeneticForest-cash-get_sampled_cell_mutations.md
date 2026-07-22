# Getting the sampled cells' mutations

This method returns the mutations of the sample cells.

## Arguments

- with_germline:

  A Boolean flag to report germline mutations too (default: `FALSE`).

## Value

A data frame reporting `cell_id`, `chr`, (i.e., the mutation
chromosome), `from` (i.e., position in the chromosome), `allele` (in
which the mutation occurs), `ref`, `alt`, `type` (i.e., either `"SNV"`
or `"indel"`), `cause`, and `class` (i.e., `"driver"`, `"passenger"`,
`"germinal"` or `"pre-neoplastic"`) for each mutation in the sampled
cell genomes.

## Details

This method builds a data frame representing all the SNV and the indel
mutations in the cells sampled during the simulation and represented by
the leaves of the phylogenetic forest. The data frame also reports the
allele in which the mutations occur to support double occurrences due to
CNAs.

## See also

[`PhylogeneticForest$get_sampled_cell_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_CNAs.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the sampled cell mutations
mutations <- forest$get_sampled_cell_mutations()

# print the first lines of the data frame
head(mutations)
#>   chr     from allele            ref  alt cause         nature cell_id
#> 1  22 16085675      0      GCCTCCCGA    G     A         driver  183352
#> 2  22 16095655      0              T    A  SBS1 pre-neoplastic  183352
#> 3  22 16099091      0              G GATA   ID1 pre-neoplastic  183352
#> 4  22 16165397      0              T  TTG   ID1 pre-neoplastic  183352
#> 5  22 16241450      0              T    G  SBS1 pre-neoplastic  183352
#> 6  22 16288596      0 CGCGTGCGGCGTGC    C   ID1 pre-neoplastic  183352
```
