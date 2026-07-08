# An error-free Illumina sequencer class

This class implements a perfect Illumina sequencers that does not commit
errors.

This method builds an error-free Illumina sequencer model.

## Value

A new error-free Illumina sequencer.

## See also

[`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md),
[`simulate_normal_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_normal_seq.md),
and
[`vignette("sequencing")`](https://caravagnalab.github.io/ProCESS/1.3/articles/sequencing.md)
for usage examples.

## Examples

``` r
# build a sequencer model
sequencer <- ErrorlessIlluminaSequencer()
sequencer
#> Error-less Illumina Sequencer (platform: "ILLUMINA")
```
