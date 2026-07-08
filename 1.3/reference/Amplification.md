# Creating a CNA amplification

This function creates a CNA amplification.

## Usage

``` r
Amplification(chr, chr_pos, len, allele = NULL, src_allele = NULL)
```

## Arguments

- chr:

  The name of the chromosome in which the CNA occurs.

- from:

  The position in the chromosome where the CNA occurs.

- len:

  The CNA length.

- allele:

  The allele in which the amplification is placed. (optional)

- src_allele:

  The allele from which the region is amplified. (optional)

## See also

[`Deletion()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Deletion.md)
to build a deletion;
[`CNA()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md)
to build both amplifications and deletions.

## Examples

``` r
# create an amplification CNA
cna <- Amplification("X", 20002, 100)

cna
#> CNA(type: "A", chr: "X", from: 20002, length: 100)
```
