# Getting the mutation cause

This method returns the mutation cause.

## Value

The mutation cause.

## Details

Every mutation is associated to a cause depending on whether it is part
of a genomic characterisation of a mutant or it is caused by a specific
profile. This method returns such a cause whenever it is available.

## See also

[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
# let us build a SNV without specifying any cause for it
snv <- SNV("X", 20002, "T", "A")

# get the cause of `snv` (i.e., "NA")
snv$get_cause()
#> [1] NA

# we can also build a SNV, specifying a cause for it
snv <- SNV("X", 20002, "T", "A", cause = "SBS13")

# get the cause of `snv` (i.e., "SBS13")
snv$get_cause()
#> [1] "SBS13"
```
