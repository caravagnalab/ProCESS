# Creating a CNA deletion

This function creates a CNA deletion.

## Usage

``` r
Deletion(chr, chr_pos, len, allele = NULL)
```

## Arguments

- chr:

  The name of the chromosome in which the CNA occurs.

- from:

  The position in the chromosome where the CNA occurs.

- len:

  The CNA length.

- allele:

  The allele in which the deletion occurs. (optional)

## See also

[`Amplification()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Amplification.md),
[`CNA()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md),
[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA_class.md)

## Examples

``` r
# create a deletion CNA
cna <- Deletion("Y", 40020, 200)

cna
#> CNA(type: "D", chr: "Y", from: 40020, length: 200)
```
