# Iterating over Forest Nodes

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS provides a practical interface to iterate over the nodes of any
forest by using the natures
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)
and
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md).

Let us reload the sample forest produced in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/sampling.md).

``` r

library("ProCESS")
#> 
#> Attaching package: 'ProCESS'
#> The following object is masked from 'package:utils':
#> 
#>     example

# load the phylogenetic forest produced by sampling.md
sample_forest <- load_forest("sample_forest.sff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# show some information about the forest
sample_forest
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 21182
#>   # of leaves: 5425
#>   samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}
```

The function
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
accepts as parameter a forest and returns either a
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)
or a
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md)
object depending on the forest type. These natures implement iterators
over the forest nodes and, at each instant, point to one of the forest
nodes. Their objects provide the properties `node` and `done`, as well
as the method [`step()`](https://rdrr.io/r/stats/step.html). The
property `node` maintains the node pointed by the iterator and is either
a
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)
or a
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md)
object depending on the type of the forest passed to
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md).
Instead, the tour property `done` is a Boolean value that is `TRUE` if
and only if the end of the iteration has been reached. Finally, the
method [`step()`](https://rdrr.io/r/stats/step.html) lets the iterator
move to the next node when available.

``` r

node_tour <- get_node_tour(sample_forest)
node_tour
#> SampleForestNodeTour(with_labels = FALSE)

# get the first node in the tour
node_tour$node
#> SampleForestNode(cell_id = 1, species = "A[E1]")

# test whether the tour has been completed
node_tour$done
#> [1] FALSE

# go to the next node
node_tour$step()

# the second value of the labelling tour
node_tour$node
#> SampleForestNode(cell_id = 2, species = "A[E1]")
```

Since the forest type is
[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
returns a
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)
object and the property `node` has type
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md).

Let us reload the phylogenetic forest produced in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)
and build a new tour.

``` r

# load the phylogenetic forest produced by mutations.md
phylo_forest <- load_forest("phylo_forest.pff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# show some information about the forest
phylo_forest
#> PhylogeneticForest
#>   # of trees: 1
#>   # of nodes: 21182
#>   # of leaves: 5425
#>   samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}
#> 
#>   # of emerged SNVs and indels: 16731
#>   # of emerged CNAs: 37

# build a tour over `phylo_forest` nodes
node_tour <- get_node_tour(phylo_forest)
node_tour
#> PhylogeneticForestNodeTour(with_labels = FALSE, with_genomes = FALSE)

# get the first node in the tour
node_tour$node
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")
```

The tour visits nodes in depth first search order. Thus, the first node
to be visited is a forest root.

``` r

# check whether the node is a forest root
node_tour$node$is_root
#> [1] TRUE

# check whether the node is a forest leaf
node_tour$node$is_leaf
#> [1] FALSE

# get the name of the sample that collected the cell corresponding to
# the node
node_tour$node$sample_name
#> [1] NA
```

#### Visiting the Forest Leaves

By setting the optional Boolean parameter `only_leaves` to `TRUE`,
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
produces an iterator over the forest leaves.

``` r

# build a tour over `phylo_forest` node
node_tour <- get_node_tour(phylo_forest, only_leaves = TRUE)

# check whether the first node is a forest leaf
node_tour$node$is_leaf
#> [1] TRUE

# move to the next leaf
node_tour$step()

# check whether the second node is a forest leaf
node_tour$node$is_leaf
#> [1] TRUE

# get the name of the sample that collected the cell corresponding to
# the node
node_tour$node$sample_name
#> [1] "S_1_1"
```

#### Labelling Nodes

The function
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
accepts as an optional parameter `labelling_functor` and generates the
node labels by applying the labelling functor to the forest nodes in a
depth-first visit order.

The labelling functor assigns a label to each node within the forest and
must have type `function(parent_label, node)` where the type of `node`
is either a
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
object or a
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
object depending on the forest type.

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
  list(cell_id = node$cell_id, species = node$species_name)
}

# a function that label each node by the list containing its cell
# identifier and the identifier of its ancestors
ancestor_species_labelling <- function(parent_label, node) {
  # when the node a root its parent label is NULL
  if (is.null(parent_label)) {
    list(node$cell_id)
  } else {
    unique(append(parent_label, node$species_name))
  }
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

When the optional parameter `labelling_functor` is provided, the label
of the pointed forest node is available in the tour additional property
`label`.

``` r

# generate a node tour labelling each node by the identifiers
# of the corresponding cells
node_tour <- get_node_tour(phylo_forest, id_labelling)

# the first node of the tour
node_tour$node
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")

# the label of the first node
node_tour$label
#> [1] 1

# go to the next node
node_tour$step()

# the label of the second node
node_tour$label
#> [1] 2

# generate a node tour labelling each node by the ids of its ancestors
node_tour <- get_node_tour(phylo_forest, ancestor_species_labelling)

# go to the first node whose associated cell belongs to "B[E1]"
while (node_tour$node$species_name != "B[E1]") {
  node_tour$step()
}

node_tour$node
#> PhylogeneticForestNode(cell_id = 207882, species = "B[E1]")

# show the first five ancestors of the node
node_tour$label[1:5]
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] "A[E1]"
#> 
#> [[3]]
#> [1] "A[E2]"
#> 
#> [[4]]
#> [1] "B[E2]"
#> 
#> [[5]]
#> [1] "B[E1]"
```

#### Producing the Simulated Genomes

Behind the labels produced by the labelling function,
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md)
can also generate the genomes of the forest nodes. This is achieved by
using the
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)’s
optional Boolean parameter `with_genomes`. When `with_genomes` is set to
`TRUE`, the tour additional property `genome` contains a
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
object representing the genome of the cell modelled by the pointed node.

``` r

# generate a node tour with genomes
node_tour <- get_node_tour(phylo_forest, with_genomes = TRUE)

node_tour$genome
#> GenomeMutations: 1 chrs 6 alleles
```

Notice that the genome of each node can by generated by using the
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
method
[`get_genome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md).
However, the asymptotic time cost of this method is linear in the number
of ancestors of the node, while the overall cost of a tour is linear in
the number of nodes. Hence, iterating over the the node/leaf genomes is
much faster than iterating over the nodes and generating their genomes
by using
[`get_genome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md).

##### Accessing Cell Genome Data

The class
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
provides methods to access cell genome data.

The method
[`GenomeMutations$get_allele_fragments()`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations-cash-get_allele_fragments.md)
returns a data frame containing all the allele fragments in cell genome.

``` r

# get the genome of the pointed node
genome <- node_tour$genome

# see the genome structure (i.e., chromosomes, alleles, allele fragments).
# Here we exclusively see chromosome 22 because of the mutation engine
# demo setup.
genome$get_allele_fragments()
#>   chr allele src.allele     from     size
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
#>   chr     from allele ref alt cause         nature
#> 1  22 16133982      0   T   A  SBS1 pre-neoplastic
#> 2  22 16270440      0   T   A  SBS1 pre-neoplastic
#> 3  22 16277531      0   C  CT   ID1 pre-neoplastic
#> 4  22 16391790      0   A   C  SBS1 pre-neoplastic
#> 5  22 16440524      0   T   C  SBS1 pre-neoplastic
#> 6  22 16507752      0   A   T  SBS1 pre-neoplastic
```

Users can also extract a fragment of the cell genome. For instance, the
following snippet generates a
[`GenomeFragment`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeFragment.md)
object representing the fragment beginning at position 16085655,
chromosome 22, and allele 1 whose length is 100 nucleotides.

``` r

# extract a fragment from the genome
fragment <- genome$get_fragment("22", 1, 16085655, 100)

fragment
#> chr22(1)[16085655-16085755]
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
#> [1] 1
#> 
#> $from
#> [1] 16085655
#> 
#> $size
#> [1] 108

# the fragment DNA sequence
fragment$sequence
#> [1] "CATGCCATTCTCCTGCCTCAGGTAGCTGGGACTACTGGCACCCACCACTACGCCCGGATACTTTTTGTATTTTTAGTAGAGACGGGGTTTCACTGTGTTA"

# the mutations laying in the fragment
fragment$get_mutations()
#>   chr allele     from       ref alt cause nature
#> 1  22      1 16085675 GCCTCCCGA   G     A driver

# the CIGAR code for the fragment
fragment$get_CIGAR()
#> [1] "21M8D79M"
```

Li, Heng, Bob Handsaker, Alec Wysoker, et al. 2009. “The Sequence
Alignment/Map Format and SAMtools.” *Bioinformatics* 25 (16): 2078–79.
<https://doi.org/10.1093/bioinformatics/btp352>.
