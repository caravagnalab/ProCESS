# The tissue rectangle upper corner

This is the vector of the maxima among all the rectangle dimensions.

## Examples

``` r
rect <- new(TissueRectangle, c(500, 500), c(550, 550))

# get the rectangle upper corner, i.e., (550, 550)
rect$upper_corner
#> [1] 550 550
```
