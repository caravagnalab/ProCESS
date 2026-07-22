# Getting the supported setups

This method returns the supported codes for predefined set-up.

## Usage

``` r
get_mutation_engine_codes()
```

## Value

A data frame reporting the code and a description for each supported
predefined set-up.

## See also

[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
to build a mutation engine

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)

## Examples

``` r
# get the list of supported mutation engine set-up codes
get_mutation_engine_codes()
#>     name            description
#> 1 GRCh37  Homo sapiens (GRCh37)
#> 2 GRCh38  Homo sapiens (GRCh38)
#> 3   demo A demonstrative set-up
```
