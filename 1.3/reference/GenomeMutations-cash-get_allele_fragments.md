# Getting genome allele fragments

This method returns a data frame representing the allele fragments in
the genome.

## Value

A data frame consisting of 5 columns: `chr`, `allele`, `src_allele`,
`from`, and `size`. Each row represent a allele fragment. The columns
`chr`, and `allele` represent the fragment's chromosome and allele,
respectively. The column `allele_src` stores the allele from which the
allele of the fragment is derived. The columns `from` and `size`
maintain the first position of the fragment in the wild-type allele and
its size.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the genome of the cell having 2 as the identifier
genome <- forest$get_node(2)$get_genome()

# get the genome allele fragments
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
