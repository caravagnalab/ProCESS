# Getting information about the fragment aligning to the reference

This method returns information about the fragment of the current genome
that align with the provided region in the reference.

## Arguments

- chromosome_name:

  The name of the chromosome of the reference region.

- allele:

  The allele of the reference region.

- from:

  The first position in the chromosome of the reference region.

- size:

  The size of the reference region.

## Value

A named list whose names are `chr`, `allele`, `from`, and `length`
representing the fragment in the current genome that aligns on the
region in chromosome `chromosome_name`, allele `allele` from position
`from` and whose size is `size`.

## See also

[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get the region in chromosome 22 allele 0 aligning to the
# reference region from position 16085625 whose size is 100
genome$get_region_aligned_to_ref("22", 0, 16085625, 100)
#> $chr
#> [1] "22"
#> 
#> $allele
#> [1] 0
#> 
#> $from
#> [1] 16085625
#> 
#> $length
#> [1] 92
#> 
```
