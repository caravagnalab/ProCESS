# Build a new rectangle of tissue

Build a new rectangle of tissue

## Examples

``` r
# build the rectangle [500,550]x[450,475]
rect <- new(TissueRectangle, c(500, 450), c(550, 475))

rect
#> TissueRectangle((500,450),(550,475))

# build the rectangle [500,550]x[450,475]
rect <- new(TissueRectangle, c(500, 450), 50, 25)

rect
#> TissueRectangle((500,450),(549,474))
```
