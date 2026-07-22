# Representing a genome fragment

This class represents genome fragment.

## Details

This class represents genome fragment. The objects of this class are
built by
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md).
They provides the following properties and methods:

- [`get_CIGAR()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment-cash-get_CIGAR.md)
  returns the CIGAR code of the fragment with respect to the reference
  genome.

- [`get_covered_reference_region()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment-cash-get_covered_reference_region.md)
  returns the reference genome region covered by the fragment.

- [`get_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment-cash-get_mutations.md)
  returns the mutations occurring on the fragment.

- [`sequence`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment-cash-sequence.md)
  stores the fragment DNA sequence.

## See also

[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md),
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
