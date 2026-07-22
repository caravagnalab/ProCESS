# Getting the cell in which a mutation emerged

This method returns the identifier of the cell in which a mutation
occurs for the first time.

## Arguments

- mutation:

  A mutation being a SNV, a indel, or a CNA.

## Value

The identifier of the cell in which a mutation occurs for the first
time.

## See also

[`vignette("mutations")`](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# consider a mutation in the forest
mutation <- Mutation("22", 35396109, ref = "CTGA", alt = "C")
mutation
#> indel(chr: 22, from: 35396109, allele: random, ref: CTGA, alt: C)

# get the identifiers of the cells in which the mutation arose, i.e.,
# they have the mutation, but their parent have not
cell_ids <- forest$get_first_occurrences(mutation)
cell_ids
#> [[1]]
#> [1] 229
#> 

# get the corresponding node
node <- forest$get_node(cell_ids[[1]])

# the mutation is among the arising mutation in the node
node$arising_mutations
#>   order type CNA_type chr    start      end  ref alt allele src.allele
#> 1     1  SID     <NA>  22 35396109 35396112 CTGA   C      3         NA
#>      nature cause code
#> 1 passenger   ID2 <NA>
```
