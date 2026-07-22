# Labelling forest nodes

This method generates a
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)

## Usage

``` r
get_node_tour(forest, labelling_functor, init_value, only_leaves,
                      with_genomes)
```

## Arguments

- forest:

  Either a
  [`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md)
  or a
  [`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
  object.

- labelling_functor:

  Depending on the type of `forest`, a function of the type
  `label_type (*)(label_type, SampleForestNode)`, when `forest` is a
  [`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md)
  object, or `label_type (*)(label_type, PhylogeneticForestNode)`, when
  `forest` is is a
  [`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
  object (default: `NULL`).

- init_value:

  The initial value of the labelling process (default: `NULL`).

- only_leaves:

  A Boolean flag to iterate exclusively over the forest leaves (default:
  `FALSE`).

- with_genomes:

  A Boolean flag to also generate the corresponding cell genomes. It can
  be set to TRUE exclusively when `forest` is a
  [`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
  object. (default: `FALSE`)

## Value

Either a
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md)
or a
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md)
iterating over `forest`'s nodes. When `only_leaves` is set to `TRUE`,
the returned tour iterates over `forest`'s leaves.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get a tour of the forest nodes
node_tour <- get_node_tour(forest)

# the first node in the tour
node_tour$node
#> SampleForestNode(cell_id = 1, species = "A[E1]")

# the first node is not a leaf
node_tour$node$is_leaf
#> [1] FALSE

# move to the next node
node_tour$step()

# not even the second node is a leaf
node_tour$node$is_leaf
#> [1] FALSE

# use a sample forest example
forest <- example("SampleForest")

# get a tour of the forest leaves
node_tour <- get_node_tour(forest, only_leaves = TRUE)

# the first node in the tour is a leaf
node_tour$node
#> SampleForestNode(cell_id = 183352, species = "A[E2]")
node_tour$node$is_leaf
#> [1] TRUE

# move to the next node
node_tour$step()

# the second node too
node_tour$node$is_leaf
#> [1] TRUE

# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get a tour of the forest nodes with their genomes
node_tour <- get_node_tour(forest, with_genomes = TRUE)

# the first node in the tour is a leaf
node_tour$node
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")

# the forest was build using the setup "demo". Thus, the
# node genomes have only chromosome 22
node_tour$genome
#> GenomeMutations: 1 chrs 6 alleles

# use a sample forest example
forest <- example("SampleForest")

# we define a function to collect the tour labels
collect_labels <- function(tour) {
  total <- NULL

  # `tour$done` is `TRUE` iff the tour ended
  while (!tour$done) {
    if (is.null(total)) {
      # `tour$label` is the node label
      total <- list(tour$label)
    } else {
      total <- append(total, tour$label)
    }

    # `tour$step()` advances to the next node
    # in the tour
    tour$step()
  }

  total
}

print("Functor 1")
#> [1] "Functor 1"
# a labelling functor
labelling_functor1 <- function(label, node) {
  # the nodes are labelled by the identifiers of the associated cells
  node$cell_id
}

# since `labelling_functor1` does not use `label`, we can omit the
# parameter `init_value`
tour <- get_node_tour(forest, labelling_functor1)

print("Functor 1 - All nodes")
#> [1] "Functor 1 - All nodes"
print(collect_labels(tour)[1:5])
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 6
#> 
#> [[4]]
#> [1] 8
#> 
#> [[5]]
#> [1] 10
#> 

# since `labelling_functor1` does not use `label`, we can omit the
# parameter `init_value`
tour <- get_node_tour(forest, labelling_functor1, only_leaves = TRUE)

print("Functor 1 - Only leaves")
#> [1] "Functor 1 - Only leaves"
print(collect_labels(tour)[1:5])
#> [[1]]
#> [1] 183352
#> 
#> [[2]]
#> [1] 165816
#> 
#> [[3]]
#> [1] 177858
#> 
#> [[4]]
#> [1] 177859
#> 
#> [[5]]
#> [1] 197872
#> 

labelling_functor2 <- function(label, node) {
  # the nodes are labelled by their visiting order
  label + 1
}

# `labelling_functor2` uses `label` and we must specify the
# parameter `init_value`
tour <- get_node_tour(forest, labelling_functor2,
                       init_value = 0, only_leaves = TRUE)

print("Functor 2 - Only leaves")
#> [1] "Functor 2 - Only leaves"
print(collect_labels(tour)[1:5])
#> [[1]]
#> [1] 24
#> 
#> [[2]]
#> [1] 39
#> 
#> [[3]]
#> [1] 40
#> 
#> [[4]]
#> [1] 40
#> 
#> [[5]]
#> [1] 34
#> 

a <- 3.14

# this functor uses a global variable to compute the
# labels
labelling_functor3 <- function(label, node) {
  # the nodes are labelled by their visiting order multiplied
  # by the value in the global variable `a`
  label + a
}

tour <- get_node_tour(forest, labelling_functor3,
                       init_value = 0, only_leaves = TRUE)

print("Functor 3 - Only leaves")
#> [1] "Functor 3 - Only leaves"
print(collect_labels(tour)[1:5])
#> [[1]]
#> [1] 75.36
#> 
#> [[2]]
#> [1] 122.46
#> 
#> [[3]]
#> [1] 125.6
#> 
#> [[4]]
#> [1] 125.6
#> 
#> [[5]]
#> [1] 106.76
#> 

set.seed(0)

# this functor uses a random function
labelling_functor4 <- function(label, node) {
  # the nodes are randomly labelled
  label + sample(-100:100, 1)
}

tour <- get_node_tour(forest, labelling_functor4,
                       init_value = 0, only_leaves = TRUE)

print("Functor 4 - Only leaves")
#> [1] "Functor 4 - Only leaves"
print(collect_labels(tour)[1:5])
#> [[1]]
#> [1] -62
#> 
#> [[2]]
#> [1] 49
#> 
#> [[3]]
#> [1] -199
#> 
#> [[4]]
#> [1] -198
#> 
#> [[5]]
#> [1] -89
#> 
```
