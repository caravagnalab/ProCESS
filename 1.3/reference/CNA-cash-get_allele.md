# Getting the CNA allele

This method returns the identifier of the allele in which the CNA
occurred.

## Value

The allele in which CNA occurred.

## Details

If the CNA is an amplification corresponds to the new allele identifier.
If, instead, the CNA is a deletion is the identifier of the allele on
which the deletion occurred.

## See also

[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA_class.md)

## Examples

``` r
# create an amplification CNA
cna <- Amplification("X", 20002, 100, 1, 0)

# get the allele in which `cna` occurs (i.e., 1)
cna$get_allele()
#> [1] 1
```
