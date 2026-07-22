# Getting the error rate

This method returns the sequencing error rate of the simulated illumina
sequencer.

## Value

The sequencing error rate of the simulated sequencer.

## See also

[`BasicIlluminaSequencer`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer_class.md)

## Examples

``` r
# build a basic Illumina sequencer model whose errors occur
# at rate 4e-3
sequencer <- BasicIlluminaSequencer(4e-3)

sequencer$error_rate
#> [1] 0.004

sequencer$error_rate <- 5e-2

sequencer$error_rate
#> [1] 0.05
```
