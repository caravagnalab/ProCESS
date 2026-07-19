# Get the available data structure examples

This function returns the available data structure examples.

## Usage

``` r
available_examples()
```

## Value

A data frame describing the available examples and consisting of two
columns: the names of the data structure examples (column `name`) and
their descriptions (column `description`).

## See also

[`example()`](https://caravagnalab.github.io/ProCESS/1.3/reference/example.md)

## Examples

``` r
# get the data frame of the available data structure examples
available_examples()
#>                                name
#> 1                      SampleForest
#> 2                PhylogeneticForest
#> 3       SampleForest - no epistates
#> 4 PhylogeneticForest - no epistates
#> 5                Sequencing results
#>                                                                          description
#> 1       A sample forest with two mutations, two epigenetic states, and four samples.
#> 2 A phylogenetic forest with two mutations, two epigenetic states, and four samples.
#> 3        A sample forest with two mutations, no epigenetic states, and four samples.
#> 4  A phylogenetic forest with two mutations, no epigenetic states, and four samples.
#> 5                The result of a 10x sequencing of the example "PhylogeneticForest".
```
