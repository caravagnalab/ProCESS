# Getting the genome information

This method returns information about the genome.

## Details

This method returns a data frame reporting the name (column `name`), the
size (column `size`), and the number of alleles (column
`num_of_alleles`) of each chromosome.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# get the genome information
m_engine$get_genome_info()
#>   chr     size num_of_alleles
#> 1  22 51304566              2
```
