# Building a basic Illumina sequencer simulator

This method builds a basic Illumina sequencer model.

## Usage

``` r
BasicIlluminaSequencer(error_rate, random_quality_scores)
```

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

## See also

[`BasicIlluminaSequencer`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer_class.md)

## Examples

``` r
# build a sequencer model having error rate 4e-3
sequencer <- BasicIlluminaSequencer(error_rate=4e-3)
sequencer
#> Basic Illumina Sequencer (platform: "ILLUMINA" error rate: 0.004000 random quality scores)
```
