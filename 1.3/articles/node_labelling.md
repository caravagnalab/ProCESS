# Forest Node Labelling

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS supports recursive labelling of forest nodes. The function
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md)
accepts a forest and a labelling function as input and generates an
iterator over the forest node labels by applying the labelling function
to the forest nodes in a depth-first visit order. The labelling function
takes as input the label of the parent node, if available, and a node,
and computes the label of the node itself.

#### Labelling Functions

The labelling function assigns a label to each node within the forest
and must have type `function(parent_label, node)` where the type of
`node` is either a
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
object or a
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
object depending on the forest type. These two classes expose properties
of cells represented by forest nodes. For instance, the cell identifiers
or the cell species. The following examples only report some of the
available properties. Please, refer to
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
and
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
for a complete overview.

A labelling function can label each node with the identifier of the
corresponding cell.

``` r

# a function that label each node by the identifier of the corresponding cell
id_labelling <- function(label, node) {
  node$cell_id
}

# a function that label each node by a list of the identifier and species
# of the corresponding cell
id_species_labelling <- function(parent_label, node) {
  list(cell_id = node$cell_id, species = node$species)
}

# a function that label each node by the list containing its cell
# identifier and the identifier of its ancestors
ancestor_id_labelling <- function(parent_label, node) {
  # when the node a root its parent label is NULL
  if (is.null(parent_label)) {
    list(node$cell_id)
  } else {
    append(parent_label, node$cell_id)
  }
}

# a function to label each node by a list consisting of two Boolean values. The
# first is named `leaf` and is `TRUE` if and only if the node is a leaf of the
# forest. The second is named `root` and is `TRUE` if and only if the node is a
# root.
node_prop_labelling <- function(parent_label, node) {
  list(leaf = node$is_leaf, root = node$is_root)
}
```

The labelling function can also access to external variables. For
instance, we can label forest nodes by using values sampled in a list by
using the following labelling function.

``` r

values <- 1:100

# this labelling function labels each node by randomly sampling the
# content of values 
random_labelling <- function(parent_label, node) {
  sample(values, 1)
}
```

#### Iterating over Forest Nodes

Let us reload the phylogenetic forest produced in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md).

``` r

library(ProCESS)

# load the phylogenetic forest produced by mutations.md
phylo_forest <- load_phylogenetic_forest("phylo_forest.sff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# show some information about the forest
phylo_forest
#> PhylogeneticForest
#>   # of trees: 1
#>   # of nodes: 21182
#>   # of leaves: 5425
#>   samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}
#> 
#>   # of emerged SNVs and indels: 16669
#>   # of emerged CNAs: 54
```

Each node of the forest represents a cell that can be labelled by the
function
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md).
This function returns either a
[`SampleForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestLabelTour.md)
or a
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
object depending on the type of the forest parameter. These classes
implement iterators over the forest nodes and, at each instant, refer to
one of the nodes. They provide the properties `value` and `done`, as
well as the method [`step()`](https://rdrr.io/r/stats/step.html). The
property `value` is a named list consisting of both the identifier of
the cell corresponding to the currently referred node (name `cell_id`)
and the node label (name `label`). Instead, the property `done` is a
Boolean value that is `TRUE` if and only if the end of the iteration has
been reached. Finally, the method
[`step()`](https://rdrr.io/r/stats/step.html) lets the iterator move to
the next node when available.

``` r

# generate a label tour labelling each node by the identifiers
# of the corresponding cells
label_tour <- get_label_tour(phylo_forest, id_labelling)

# the first value of the labelling tour
label_tour$value
#> $cell_id
#> [1] 1
#> 
#> $label
#> [1] 1

# test whether the tour has been completed
label_tour$done
#> [1] FALSE

# go to the next node
label_tour$step()

# the second value of the labelling tour
label_tour$value
#> $cell_id
#> [1] 2
#> 
#> $label
#> [1] 2

# generate a label tour labelling each node by its properties
label_tour <- get_label_tour(phylo_forest, node_prop_labelling)

# the first value of the labelling tour. Notice that the first
# node is a root, not a leaf
label_tour$value
#> $cell_id
#> [1] 1
#> 
#> $label
#> $label$leaf
#> [1] FALSE
#> 
#> $label$root
#> [1] TRUE
```

#### Iterating over Forest Leaves

The function
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md)
has the optional Boolean parameter `only_leaves` that, when set to
`TRUE` allows to iterate over the leaves.

``` r

# generate a label tour labelling each forest leaf by its properties
label_tour <- get_label_tour(phylo_forest, node_prop_labelling,
                             only_leaves = TRUE)

# since the parameter `only_leaves` is set to `TRUE`, in this case, the first
# node is not a root, but a leaf
label_tour$value
#> $cell_id
#> [1] 23281
#> 
#> $label
#> $label$leaf
#> [1] TRUE
#> 
#> $label$root
#> [1] FALSE

# ... so are all labels in the tour
label_tour$step()

label_tour$value
#> $cell_id
#> [1] 56990
#> 
#> $label
#> $label$leaf
#> [1] TRUE
#> 
#> $label$root
#> [1] FALSE
```

### Producing the Simulated Genomes

Behind the labels produced by the labelling function,
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
can also label each node by the genome of the cell corresponding to the
referred node. This is achieved by using the optional Boolean parameter
`with_genome` of the function
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md).
When `with_genome` is set to `TRUE`, the property `value` contains a
further element, whose name is `genome`, representing the cell’s genome.
In particular, this value is a
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
object.

``` r

# generate a label tour labelling each forest leaf by its properties
label_tour <- get_label_tour(phylo_forest, node_prop_labelling,
                             only_leaves = TRUE, with_genomes = TRUE)

label_tour$value
#> $cell_id
#> [1] 23281
#> 
#> $label
#> $label$leaf
#> [1] TRUE
#> 
#> $label$root
#> [1] FALSE
#> 
#> 
#> $genome
#> GenomeMutations: 1 chrs 6 alleles
```

#### Accessing Cell Genome Data

The class
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
provides methods to access cell genome data.

The method
[`GenomeMutations$get_allele_fragments()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_allele_fragments.md)
returns a data frame containing all the allele fragments in cell genome.

``` r

# get the genome of the current label tour value
genome <- label_tour$value$genome

# see the genome structure (i.e., chromosomes, alleles, allele fragments).
# Here we exclusively see chromosome 22 because of the mutation engine
# demo setup.
genome$get_allele_fragments()
#>   chr allele src_allele     from     size
#> 1  22      0          0        1 51304566
#> 2  22      1          1        1  5009999
#> 3  22      1          1  5210000 46094567
#> 4  22      2          2 10303470   200000
#> 5  22      3          0        1 51304566
#> 6  22      4          1        1  5009999
#> 7  22      4          1  5210000 46094567
#> 8  22      5          2 10303470   200000
```

The mutations accumulated in the cell genome can be obtained calling the
method
[`GenomeMutations$get_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_mutations.md).

``` r

# get the first mutations in the genome
head(genome$get_mutations())
#>   chr allele     from ref alt causes        classes
#> 1  22      0 16133982   T   A   SBS1 pre-neoplastic
#> 2  22      0 16270440   T   A   SBS1 pre-neoplastic
#> 3  22      0 16277531   C  CT    ID1 pre-neoplastic
#> 4  22      0 16391790   A   C   SBS1 pre-neoplastic
#> 5  22      0 16440524   T   C   SBS1 pre-neoplastic
#> 6  22      0 16507752   A   T   SBS1 pre-neoplastic
```

Users can also extract a fragment of the cell genome. For instance, the
following snippet generates a
[`GenomeFragment`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment.md)
object representing the fragment beginning at position 16085655,
chromosome 22, and allele 0 whose length is 100 nucleotides.

``` r

# extract a fragment from the genome
fragment <- genome$get_fragment("22", 0, 16085655, 100)

fragment
#> chr22(0)[16085655-16085755]
```

#### Genome Fragment Data

Behind its position, any
[`GenomeFragment`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment.md)
object maintains the DNA sequence corresponding to the fragment, the
mutations occurring to it, and its CIGAR code (Li et al. 2009) with
respect to the reference genome.

``` r

# the region covered in the reference genome
fragment$get_covered_reference_region()
#> $chr
#> [1] "22"
#> 
#> $allele
#> [1] 0
#> 
#> $from
#> [1] 16085655
#> 
#> $size
#> [1] 100

# the fragment DNA sequence
fragment$get_sequence()
#> [1] "CATGCCATTCTCCTGCCTCAGGTAGCTGGGACTACTGGCACCCACCACTACGCCCGGATACTTTTTGTATTTTTAGTAGAGACGGGGTTTCACTGTGTTA"

# the mutations laying in the fragment
fragment$get_mutations()
#> [1] chr     allele  from    ref     alt     causes  classes
#> 1  22      0 16085675 GCCTCCCGA   G      A  driver

# the CIGAR code for the fragment
fragment$get_CIGAR()
#> [1] "21M8D79M"
```

#### Iterating over Genomes without Labels

The function
[`get_genome_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_genome_tour.md)
takes as input a
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
object and returns an iterator that produce the simulated cell genomes
without any further node labels. The returned value is a
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
whose property `value` does not have the `label` element.

``` r

genome_tour <- get_genome_tour(phylo_forest)

genome_tour$value
#> $cell_id
#> [1] 1
#> 
#> $genome
#> GenomeMutations: 1 chrs 6 alleles

genome_tour$step()

genome_tour$value
#> $cell_id
#> [1] 2
#> 
#> $genome
#> GenomeMutations: 1 chrs 6 alleles
```

As in the case of
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md),
the function
[`get_genome_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_genome_tour.md)
as an optional Boolean parameter `only_leaves` that, when set to `TRUE`,
allows to focus on the genomes of forest leaves.

``` r

genome_tour <- get_genome_tour(phylo_forest, only_leaves = TRUE)

genome_tour$value
#> $cell_id
#> [1] 23281
#> 
#> $genome
#> GenomeMutations: 1 chrs 6 alleles

genome_tour$step()

genome_tour$value
#> $cell_id
#> [1] 56990
#> 
#> $genome
#> GenomeMutations: 1 chrs 6 alleles
```

Li, Heng, Bob Handsaker, Alec Wysoker, et al. 2009. “The Sequence
Alignment/Map Format and SAMtools.” *Bioinformatics* 25 (16): 2078–79.
<https://doi.org/10.1093/bioinformatics/btp352>.
