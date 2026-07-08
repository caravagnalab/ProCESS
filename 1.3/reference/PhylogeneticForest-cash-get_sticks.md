# Computing the forest sticks

This method computes the sticks of the forest.

## Arguments

- birth_threshold:

  The maximum birth time for the cells associated to the returned sticks
  (optional).

## Value

The list of the forest sticks whose associated cells have birth time
smaller than or equal to `birth_threshold`. Each stick is represented as
the list of cell identifiers labelling the nodes in the stick from the
higher to the deeper in the forest.

## Details

A *crucial node* of a forest is a root of the forest, a node whose
parent belongs to a different species, or the most recent common
ancestor of two crucial nodes.

A *stick* is a path of the forest in which the only crucial nodes are
the first and the last one.

This method returns the list of the forest sticks. Each stick is
represented by the sequence of cell identifiers labelling the nodes in
the stick.

## See also

[`SampleForest$get_sticks()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_sticks.md)
for usage examples.
