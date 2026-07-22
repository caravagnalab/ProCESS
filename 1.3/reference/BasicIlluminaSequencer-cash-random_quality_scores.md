# Check the non-constant quality score model

This field is set to `TRUE` if and only if the sequencers implements a
non-constant quality score model.

## Value

`TRUE` if and only if the sequencers sequencers implements a
non-constant quality score model.

## See also

[`BasicIlluminaSequencer`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer_class.md)

## Examples

``` r
# build a basic Illumina sequencer model whose quality scores are
# non-constant
sequencer <- BasicIlluminaSequencer(4e-3)

sequencer$random_quality_scores
#> [1] TRUE

sequencer$random_quality_scores <- FALSE

sequencer$random_quality_scores
#> [1] FALSE
```
