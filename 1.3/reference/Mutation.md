# Either an SBS or an indel

This function creates SNVs and indels.

## Arguments

- chr:

  The name of the chromosome in which the indel occurs.

- from:

  The position in the chromosome where the indel occurs.

- ref:

  The reference sequence.

- alt:

  The mutation altered sequence.

- allele:

  The allele in which the mutation must occur (optional).

- cause:

  The cause of the mutation (optional).

## Details

This function generalizes the function
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)
by constructing SNVs and indels. However, it necessitates the
specification of the reference sequence, whereas
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)
can infer it from the reference sequence itself.

Another distinction between this function and
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)
lies in the order of the `ref`-`alt` parameter: in
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md),
the alt parameter precedes the optional ref parameter, while
`Mutation()` adopts the reverse order.

## See also

[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)
for SNV creation.

## Examples

``` r
# create a deletion without specifying the cause
mutation <- Mutation("X", 20002, "TAC", "T")
mutation
#> indel(chr: X, from: 20002, allele: random, ref: TAC, alt: T)

# create an insertion and do not specify the cause
mutation <- Mutation("X", 20002, "A", "AT")
mutation
#> indel(chr: X, from: 20002, allele: random, ref: A, alt: AT)

# create an insertion that must be place in allele 1
mutation <- Mutation("X", 20002, "A", "AT", allele = 1)
mutation
#> indel(chr: X, from: 20002, allele: 1, ref: A, alt: AT)

# create an insertion with a cause
mutation <- Mutation("X", 20002, "A", "AT", cause = "SBS1")
mutation
#> indel(chr: X, from: 20002, allele: random, ref: A, alt: AT, cause: "SBS1")
```
