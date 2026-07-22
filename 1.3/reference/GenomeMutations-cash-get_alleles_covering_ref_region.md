# Getting the alleles covering a reference region

This method returns the identifiers of the alleles that containing the
specified region of the reference genome

## Arguments

- chromosome_name:

  The name of the chromosome of the reference region.

- from:

  The first position in the chromosome of the reference region.

- size:

  The size of the reference region.

## Value

A list of allele identifiers. Each identifier in the list corresponds to
an allele containing the specified reference region.

## See also

[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# the genome has 6 alleles
genome
#> GenomeMutations: 1 chrs 6 alleles

# get the alleles in the genome which covers the specified region
genome$get_alleles_covering_ref_region("22", 16085625, 100)
#> [1] 0 1 3 4
```
