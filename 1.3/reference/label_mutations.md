# Label mutations using phylogenetic forest data

The function labels mutations using data about the cell in which it
occurs for the first time.

## Usage

``` r
label_mutations(seq_results, phylo_forest)
```

## Arguments

- seq_results:

  The output of
  [`simulate_seq()`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md).

- phylo_forest:

  The phylogenetic forest from which the sequencing was simulated.

## Value

A copy of the data frame `seq_results` added with the identifier of the
cell in which the mutation occurs for the first time (column
"`cell_id`"), the identifier of its ancestor (column "`ancestor`"), its
mutant (column "`mutant`"), its epigenetic state (column "`epistate`"),
its birth time (column "`birth_time`"), the sample that collected the
cell whenever available (column "`sample`"), and a cell classification
based on phylogenetic sticks (column "`label`").

## Examples

``` r
# use the phylogenetic forest example and its 10x sequencing results
forest <- example("PhylogeneticForest")
seq_results <- example("Sequencing results")

# label filtered mutations using phylogenetic forest data
labels <- label_mutations(seq_results$mutations, forest)
labels
#> # A tibble: 2,365 × 29
#>    chr       from ref    alt   cause nature S_1_1.NV S_1_1.DP S_1_1.VAF S_1_2.NV
#>    <chr>    <int> <chr>  <chr> <chr> <chr>     <int>    <int>     <dbl>    <int>
#>  1 22    16066874 G      T     SBS1  pre-n…        8       10     0.8          5
#>  2 22    16080730 GACT   G     ID1   pre-n…        4        9     0.444        5
#>  3 22    16085675 GCCTC… G     A     driver        4        8     0.5          7
#>  4 22    16095655 T      A     SBS1  pre-n…        5       11     0.455        2
#>  5 22    16099091 G      GATA  ID1   pre-n…        5       13     0.385        7
#>  6 22    16099402 A      C     SBS1  pre-n…        6        9     0.667        7
#>  7 22    16119735 A      G     SBS1  pre-n…        6       11     0.545        2
#>  8 22    16122890 C      A     SBS5  passe…        0       16     0            0
#>  9 22    16127180 T      G     SBS1  pre-n…        5        8     0.625        5
#> 10 22    16142082 A      T     SBS3  passe…        0       12     0            0
#> # ℹ 2,355 more rows
#> # ℹ 19 more variables: S_1_2.DP <int>, S_1_2.VAF <dbl>, S_2_1.NV <int>,
#> #   S_2_1.DP <int>, S_2_1.VAF <dbl>, S_2_2.NV <int>, S_2_2.DP <int>,
#> #   S_2_2.VAF <dbl>, normal.sample.NV <int>, normal.sample.DP <int>,
#> #   normal.sample.VAF <dbl>, cell_id <dbl>, ancestor <int>, depth <int>,
#> #   mutant <chr>, epistate <chr>, sample <chr>, birth_time <dbl>, label <chr>

# plotting histogram of the VAF with phylogenetic labels
plot_VAF_histogram(seq_results, labels = labels["label"],
                   cuts = c(0.02, 1))
```
