# Representing cell genome

This class represents genome mutations of phylogenetic forest's cells.

## Details

The objects of this class are built by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
and
[`PhylogeneticForestNode$get_genome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode-cash-get_genome.md).
The provide the following methods:

- [`get_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_CNAs.md)
  returns a data frame of the genome CNAs.

- [`get_allele_fragments()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_allele_fragments.md)
  returns genome allele fragments.

- [`get_alleles_covering_ref_region()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_alleles_covering_ref_region.md)
  returns the alleles covering a reference region.

- [`get_fragment()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_fragment.md)
  returns a fragment of the genome.

- [`get_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_mutations.md)
  returns a data frame of the genome mutations.

- [`get_region_aligned_to_ref()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_region_aligned_to_ref.md)
  returns information about the fragment aligning to the reference.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
