# A basic Illumina sequencer class

This class implements a basic model for Illumina sequencers.

## Details

It specifies a simulated sequencing error rate and the simulated
sequencing errors will occurs according to that rate.

The objects of this class are built by using the function
[`BasicIlluminaSequencer()`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer.md)
and provide the following fields:

- [`error_rate`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer-cash-error_rate.md)
  is the error rate of the sequencer.

- [`random_quality_scores`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer-cash-random_quality_scores.md)
  is a Boolean flag set to `TRUE` if and only if the sequencers
  implements a non-constant quality score model.

## See also

[`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md),
[`simulate_normal_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_normal_seq.md),
[`BasicIlluminaSequencer()`](https://caravagnalab.github.io/ProCESS/1.3/reference/BasicIlluminaSequencer.md),
and
[`vignette("sequencing")`](https://caravagnalab.github.io/ProCESS/1.3/articles/sequencing.md)
for usage examples.
