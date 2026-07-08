# Convert sequencing results from wide to long format

This function takes sequencing results in wide format and converts them
into a long format data frame. It extracts sample names from column
names, processes each sample separately, and then binds them together.
Finally, it renames and reorders columns to match the desired output
format.

## Usage

``` r
seq_to_long(seq_results)
```

## Arguments

- seq_results:

  A data frame containing sequencing results in wide format.

## Value

A data frame in long format with columns "`chr`", "`from`", "`to`",
"`ref`", "`alt`", "`NV`", "`DP`", "`VAF`", "`sample`", "`causes`", and
"`classes`".

## Examples

``` r
# Example data frame in wide format
seq_results <- data.frame(chr = c("chr1", "chr2"),
                          from = c(100, 200),
                          ref = c("A", "C"),
                          alt = c("T", "G"),
                          causes = c("SBS5", "SBS1"),
                          classes = c("germinal", "passenger"),
                          Sample.A.NV = c(10, 90),
                          Sample.A.DP = c(100, 100),
                          Sample.A.VAF = c(0.1, 0.9),
                          normal.sample.NV = c(45, 52),
                          normal.sample.DP = c(100, 100),
                          normal.sample.VAF = c(0.45, 0.52))
seq_results
#>    chr from ref alt causes   classes Sample.A.NV Sample.A.DP Sample.A.VAF
#> 1 chr1  100   A   T   SBS5  germinal          10         100          0.1
#> 2 chr2  200   C   G   SBS1 passenger          90         100          0.9
#>   normal.sample.NV normal.sample.DP normal.sample.VAF
#> 1               45              100              0.45
#> 2               52              100              0.52

# Convert to long format
seq_to_long(seq_results)
#>    chr from ref alt causes   classes NV  DP  VAF        sample  to
#> 1 chr1  100   A   T   SBS5  germinal 10 100 0.10      Sample.A 100
#> 2 chr2  200   C   G   SBS1 passenger 90 100 0.90      Sample.A 200
#> 3 chr1  100   A   T   SBS5  germinal 45 100 0.45 normal.sample 100
#> 4 chr2  200   C   G   SBS1 passenger 52 100 0.52 normal.sample 200
```
