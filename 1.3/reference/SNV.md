# Creating an SNV

This function creates SNVs.

## Usage

``` r
SNV(chr, chr_pos, alt, ref="?", allele = NULL, cause="")
```

## Arguments

- chr:

  The name of the chromosome in which the SNV occurs.

- chr_pos:

  The position in the chromosome where the SNV occurs.

- alt:

  The base after the mutation.

- ref:

  The base before the mutation (optional).

- allele:

  The allele in which the SNV must occur (optional).

- cause:

  The cause of the SNV (optional).

## See also

[`Mutation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md)
for building SNVs and indels,
[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
# create a SNV without specifying the cause and context
snv <- SNV("X", 20002, "T")
snv
#> SNV(chr: X, from: 20002, allele: random, ref: ?, alt: T)

# create a SNV and do not specify the cause
snv <- SNV("X", 20002, "T", "A")
snv
#> SNV(chr: X, from: 20002, allele: random, ref: A, alt: T)

# create a SNV that must be place in allele 1
snv <- SNV("X", 20002, "T", allele = 1)
snv
#> SNV(chr: X, from: 20002, allele: 1, ref: ?, alt: T)

# create a SNV with a cause
snv <- SNV("X", 20002, "T", cause = "SBS1")
snv
#> SNV(chr: X, from: 20002, allele: random, ref: ?, alt: T, cause: "SBS1")
```
