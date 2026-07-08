# Labelling forest nodes

This method generates a
[`SampleForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestLabelTour.md)

## Usage

``` r
get_label_tour(forest, labelling_functor, init_value, only_leaves,
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
  object.

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
[`SampleForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestLabelTour.md)
or a
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
iterates over `forest`'s nodes and applies labels from the root down to
the leaves. The labels are computed by using `labelling_functor`. For
each of `forest`'s nodes, `labelling_functor` takes as input the parent
label and the current node in the form of either a
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
or
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
object, depending on the types of `forest` and `labelling_functor` and
returns the label of the current node. The value `init_label` may be
used as parent label for all forest roots. The returned object
exclusively iterates over `forest`'s leaves if and only if `only_leaves`
is set to `TRUE`.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
[`SampleForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestLabelTour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md),
[`get_genome_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_genome_tour.md),
[`vignette("node_labelling")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_labelling.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add the mutant "A"
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation until "A" has less than 15 cells
sim$run_up_to_size("A", 15)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


sim$sample_cells("S_1_1", bottom_left = c(500, 500),
                 top_right = c(502, 502))

forest <- sim$get_sample_forest()

# we define a function to collect the tour labels
collect_labels <- function(tour) {
  total <- NULL

  # `SampleForestLabelTour$done` is `TRUE` iff the tour ended
  while (!tour$done) {
    if (is.null(total)) {
      # `SampleForestLabelTour$value` is a pair cell
      #  identifier for the current node and node label
      total <- tour$value
    } else {
      total <- rbind(total, tour$value)
    }

    # `SampleForestLabelTour$step()` advances to the next node
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
tour <- get_label_tour(forest, labelling_functor1)

print("Functor 1 - All nodes")
#> [1] "Functor 1 - All nodes"
print(collect_labels(tour))
#>       cell_id label
#> total 1       1    
#>       2       2    
#>       6       6    
#>       11      11   
#>       13      13   
#>       7       7    
#>       16      16   
#>       22      22   
#>       24      24   
#>       23      23   
#>       17      17   

# since `labelling_functor1` does not use `label`, we can omit the
# parameter `init_value`
tour <- get_label_tour(forest, labelling_functor1, only_leaves = TRUE)

print("Functor 1 - Only leaves")
#> [1] "Functor 1 - Only leaves"
print(collect_labels(tour))
#>       cell_id label
#> total 13      13   
#>       24      24   
#>       23      23   
#>       17      17   

labelling_functor2 <- function(label, node) {
  # the nodes are labelled by their visiting order
  label + 1
}

# `labelling_functor2` uses `label` and we must specify the
# parameter `init_value`
tour <- get_label_tour(forest, labelling_functor2,
                       init_value = 0, only_leaves = TRUE)

print("Functor 2 - Only leaves")
#> [1] "Functor 2 - Only leaves"
print(collect_labels(tour))
#>       cell_id label
#> total 13      5    
#>       24      6    
#>       23      5    
#>       17      4    

a <- 3.14

# this functor uses a global variable to compute the
# labels
labelling_functor3 <- function(label, node) {
  # the nodes are labelled by their visiting order multiplied
  # by the value in the global variable `a`
  label + a
}

tour <- get_label_tour(forest, labelling_functor3,
                       init_value = 0, only_leaves = TRUE)

print("Functor 3 - Only leaves")
#> [1] "Functor 3 - Only leaves"
print(collect_labels(tour))
#>       cell_id label
#> total 13      15.7 
#>       24      18.84
#>       23      15.7 
#>       17      12.56

set.seed(0)

# this functor uses a random function
labelling_functor4 <- function(label, node) {
  # the nodes are randomly labelled
  label + sample(-100:100, 1)
}

tour <- get_label_tour(forest, labelling_functor4,
                       init_value = 0, only_leaves = TRUE)

print("Functor 4 - Only leaves")
#> [1] "Functor 4 - Only leaves"
print(collect_labels(tour))
#>       cell_id label
#> total 13      39   
#>       24      64   
#>       23      110  
#>       17      -13  
```
