# Getting genome fragment

This method returns a fragment of the genome.

## Arguments

- chromosome_name:

  The name of the chromosome of the aimed fragment.

- allele:

  The allele of the aimed fragment.

- from:

  The first position in the chromosome of the aimed fragment.

- size:

  The size of the aimed fragment.

- reference_fragment:

  A reference fragment (optional).

- fragment_offset:

  The offset of the reference fragment with respect to the chromosome
  first position (optional).

## Value

The genome fragment matching the specifications.

## See also

[`GenomeFragment`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment.md)
[`GenomeMutations$get_region_aligned_to_ref()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_region_aligned_to_ref.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get a fragment
genome$get_fragment("22", 0, 16085625, 100)
#> chr22(0)[16085625-16085725]
```
