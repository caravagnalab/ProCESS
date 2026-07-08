# A basic Illumina sequencer class

This class implements a basic model for Illumina sequencers.

This method builds a basic Illumina sequencer model.

## Arguments

- error_rate:

  The error rate of the sequencer model.

- random_quality_scores:

  A Boolean flag to enable a basic non-constant quality score model.
  When it is set to `FALSE`, all the bases with no sequencing errors
  have the same quality score. The random quality score model increases
  the computation time of about 70%. (default: `TRUE`)

## Value

A basic Illumina sequencer model.

## Details

It specifies a simulated sequencing error rate and the simulated
sequencing errors will occurs according to that rate.

## See also

[`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md),
[`simulate_normal_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_normal_seq.md),
and
[`vignette("sequencing")`](https://caravagnalab.github.io/ProCESS/1.3/articles/sequencing.md)
for usage examples.

## Examples

``` r
# build a sequencer model having error rate 4e-3
sequencer <- BasicIlluminaSequencer(error_rate=4e-3)
sequencer
#> Basic Illumina Sequencer (platform: "ILLUMINA" error rate: 0.004000 random quality scores)
```
