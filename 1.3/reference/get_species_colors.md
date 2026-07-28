# Get species color maps

This function returns a species color map

## Usage

``` r
get_species_colors(object, pal_name = "Dark2", min_mutants = 4)
```

## Arguments

- object:

  The tissue simulation, the forest, or the data frame whose species
  color map is required.

- pal_name:

  The `RColorBrewer` palette name used to generate the species color map
  (default: `Dark2`).

- min_mutants:

  The minimum number of colors generated.

## Value

A named list whose name are the names of the species in `object` and
whose values are the associated colors.

## Details

This function returns a color maps for the species represented in a
tissue, in a forest, or in a data frame containing mutants and,
possible, epigenetic states.

## Examples

``` r
# get an example of `SampleForest` object
forest <- example("SampleForest")

# get the species info
forest$get_species_info()
#>   mutant epistate
#> 1      A       E1
#> 2      A       E2
#> 3      B       E1
#> 4      B       E2

# get the color map for the forest species
get_species_colors(forest)
#>     A[E1]     A[E2]     B[E1]     B[E2] 
#> "#52C098" "#1B9E77" "#FF8551" "#D95F02" 
```
