# Creating a CNA

This function creates a CNA.

## Usage

``` r
CNA(type, chr, chr_pos, len, allele = NULL, src_allele = NULL)
```

## Arguments

- type:

  The CNA type: either `"A"` or `"D"` for amplification and deletion,
  respectively.

- chr:

  The name of the chromosome in which the CNA occurs.

- from:

  The position in the chromosome where the CNA occurs.

- len:

  The CNA length.

- allele:

  The allele in which the CNA occurs. (optional)

- src_allele:

  The allele from which the region is amplified. (optional, for
  amplification only)

## See also

[`Amplification()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Amplification.md)
to build an amplification;
[`Deletion()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Deletion.md)
to build a deletion.

## Examples

``` r
# create an amplification
cna <- CNA("A", "X", 20002, 100)

cna
#> CNA(type: "A", chr: "X", from: 20002, length: 100)

# create a deletion from the allele 0
cna <- CNA("D", "Y", 101310, 205, allele = 0)

cna
#> CNA(type: "D", chr: "Y", from: 101310, length: 205, allele: 0)
```
