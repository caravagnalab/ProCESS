# The tissue rectangle lower corner

This is the vector of the minima among all the rectangle dimensions.

## Examples

``` r
# create a 51x51-rectangle from (500, 500) to (550, 550)
rect <- new(TissueRectangle, c(500, 500), c(550, 550))

# get the rectangle lower corner, i.e., (500, 500)
rect$lower_corner
#> [1] 500 500
```
