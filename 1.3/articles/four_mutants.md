# A four-mutants example

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

We build a tumour with a single mutant `A` and grow it up to 10 cells.

``` r

library("ProCESS")
#> 
#> Attaching package: 'ProCESS'
#> The following object is masked from 'package:utils':
#> 
#>     example

# set the seed of the random number generator for repeatability
set.seed(0)

sim <- TissueSimulation(width = 600, height = 600)

# set the delta time between two time series samples
sim$history_delta <- 10

# avoid drift
sim$death_activation_level <- 50

sim$add_mutant("A", c(duplication = 0.12, death = 0.05))

sim$place_cell("A", 300, 300)
sim$run_up_to_size("A", 10)
#>  [████████████████████████████████████████] 100% [00m:00s] Saving snapshot
```

We add a new mutant `B` as a sub-clone of `A` and let the simulation
evolve until `B` reaches 30 cells.

``` r

sim$add_mutant("B", c(duplication = 0.145, death = 0.06))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")

sim$run_up_to_size("B", 30)
#>  [████████████████████████████████████████] 100% [00m:00s] Saving snapshot
```

We add a further mutant `C` as a sub-clones for `A` and let the tumour
grow again up to the point at which `C` consists of 25000 cells.

``` r

sim$add_mutant("C", c(duplication = 0.15, death = 0.06))
sim$mutate_progeny(sim$choose_border_cell_in("A"), "C")

sim$run_up_to_size("C", 25000)
#>  [████████--------------------------------] 19% [00m:00s] Cells: 10994 [███████████████-------------------------] 36% [00m:00s] Cells: 19798 [███████████████████---------------------] 46% [00m:01s] Cells: 25483 [████████████████████████----------------] 58% [00m:02s] Cells: 31206 [████████████████████████████------------] 68% [00m:03s] Cells: 36976 [████████████████████████████████--------] 79% [00m:04s] Cells: 41954 [████████████████████████████████████----] 88% [00m:05s] Cells: 46679 [███████████████████████████████████████-] 95% [00m:06s] Cells: 50474 [████████████████████████████████████████] 100% [00m:07s] Saving snapshot
```

Then, we plot the tissue and the simulation state.

``` r

plot_tissue(sim)
```

![The tissue of a tumour having three mutants: A, B, and C. Both B, and
C are sub-clones of A. The three mutants have different grow rates and
appears at different simulated
time.](four_mutants_files/figure-html/unnamed-chunk-5-1.png)

``` r

plot_state(sim)
```

![The current ratio of cells per species over all tumour
cells.](four_mutants_files/figure-html/unnamed-chunk-6-1.png)

## Tissue sampling level

We can collect two samples: “S1” and “S2”.

``` r

sim$sample_cells("S1", c(145, 230), c(215, 300))
sim$sample_cells("S2", c(350, 300), c(420, 370))
```

We can plot the tissue simulation after the sampling. We label the
sampled region to improve readability, but it is not mandatory.

``` r

plot_tissue(sim) +
  ggplot2::annotate("text", x = (145 + 215) / 2, y = (230 + 300) / 2,
                    label = "S1", parse = TRUE) +
  ggplot2::annotate("text", x = (350 + 420) / 2, y = (300 + 370) / 2,
                    label = "S2", parse = TRUE)
```

![The three-mutants tumour just after collecting two samples: S1 and S2.
The former mainly contains C cells; the latter A
cells.](four_mutants_files/figure-html/unnamed-chunk-8-1.png)

After the sampling, we can add a new mutant `D` as a sub-clone of `B`,
let the simulation continue to evolve until the sum of the cardinalities
of mutants `C` and `D` is 100000, and sample the tissue again.

``` r

sim$add_mutant(name = "D", c(duplication = 0.8, death = 0.01))
sim$mutate_progeny(sim$choose_border_cell_in("B"), "D")

sim$run_until(sim$var("C") + sim$var("D") == 1e5)
#>  [██--------------------------------------] 3% [00m:00s] Cells: 49479 [████------------------------------------] 8% [00m:00s] Cells: 56707 [███████---------------------------------] 16% [00m:01s] Cells: 64651 [██████████------------------------------] 23% [00m:02s] Cells: 71525 [█████████████---------------------------] 31% [00m:03s] Cells: 79775 [█████████████████-----------------------] 40% [00m:04s] Cells: 88291 [███████████████████---------------------] 47% [00m:05s] Cells: 94379 [███████████████████████-----------------] 57% [00m:06s] Cells: 103264 [███████████████████████████-------------] 66% [00m:07s] Cells: 110599 [███████████████████████████████---------] 75% [00m:08s] Cells: 118476 [███████████████████████████████████-----] 86% [00m:09s] Cells: 127440 [███████████████████████████████████████-] 97% [00m:10s] Cells: 136800 [████████████████████████████████████████] 100% [00m:11s] Saving snapshot

sim$sample_cells("S3", c(360, 80), c(430, 150))

plot_tissue(sim) +
  ggplot2::annotate("text", x = (360 + 430) / 2, y = (80 + 150) / 2,
                    label = "S3", parse = TRUE)
```

![After collecting the samples S1 and S2, a more aggressive mutant D
arises from B and a further sample S3, containing both D and C cells, is
collected.](four_mutants_files/figure-html/unnamed-chunk-9-1.png)

``` r


plot_state(sim)
```

![The simulation state after collecting
S3.](four_mutants_files/figure-html/unnamed-chunk-9-2.png)

The time-series plot represents the species cardinalities along the
simulation.

``` r

plot_timeseries(sim)
```

![The time series plot represent the species cardinalities along the
simulation.](four_mutants_files/figure-html/unnamed-chunk-10-1.png)

The Muller plot, instead, shows the percentage of cells in each species
over the course of the simulation.

``` r

plot_muller(sim)
```

![The Muller plot overviews the percentage of cells in each species
along the simulation
evolution.](four_mutants_files/figure-html/unnamed-chunk-11-1.png)

We can build the sample forest and plot it.

``` r

sample_forest <- sim$get_sample_forest()

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

plot_forest(sample_forest) %>%
  annotate_forest(sample_forest)
```

![The sample forest of the collected
samples.](four_mutants_files/figure-html/unnamed-chunk-12-1.png)

## Cell genetic characterisation level

We first build the mutation engine to label each cell in the sample
forest. The mutation engine pre-defined setups automatically download
all the needed files. All the pre-configured setups, but “demo” requires
a COSMIC account to download the signature files from the COSMIC site.

``` r

# build the mutation engine according to the pre-defined setup "demo"
m_engine <- MutationEngine(setup_code="demo",
                           COSMIC_account = list(email = "foo@bar.com",
                                                 password = "#########"))
#>  [█---------------------------------------] 0% [00m:00s] Loading context index [████████████████████████████████████████] 100% [00m:00s] Context index loaded
#>  [█---------------------------------------] 0% [00m:00s] Loading RS index [████████████----------------------------] 28% [00m:01s] Loading RS index [███████████████████████-----------------] 55% [00m:02s] Loading RS index [█████████████████████████████-----------] 72% [00m:03s] Loading RS index [████████████████████████████████████████] 100% [00m:04s] RS index loaded
#>  [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded
```

We must genomically characterise the tissue simulation mutants,
providing passenger mutation rates and the list of driver mutations.

``` r

# genetically characterise the mutants
m_engine$add_mutant("A", passenger_rates = c(SNV = 2e-9, indel = 2e-9),
                    drivers = list(Mutation("22", 16085675, "GCCTCCCGA", "G"),
                                   CNA("D", "22", 22453799, 200000,
                                       allele = 1)))
#>  [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs [█---------------------------------------] 0% [00m:00s] Found 22 [█---------------------------------------] 0% [00m:00s] Reading 22 [█████████████---------------------------] 32% [00m:01s] Reading 22 [███████████████████████████-------------] 66% [00m:02s] Reading 22 [███████████████████████████████████████-] 97% [00m:03s] Reading 22 [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated

m_engine$add_mutant("B", passenger_rates = c(SNV = 1e-9, indel = 2e-9,
                                             CNA = 1e-9),
                    drivers = list(SNV("22", 16050185, "C", allele = 1),
                                   CNA(type = "A", chr = "22",
                                       from = 16485130, len = 200000)))
#>  [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated

m_engine$add_mutant("C", passenger_rates = list(SNV = 2e-8, indel = 2e-9),
                    drivers = list(SNV("22", 32786322, "G"),
                                   list("DGCR8 P26L", allele = 1)))
#>  [█---------------------------------------] 0% [00m:00s] Retrieving "C" SIDs [████████████████████████████████████████] 100% [00m:00s] "C"'s SIDs validated

m_engine$add_mutant("D", passenger_rates = list(SNV = 2e-9, indel = 2e-9),
                    drivers = list(SNV("22", 51240420, "T")))
#>  [█---------------------------------------] 0% [00m:00s] Retrieving "D" SIDs [████████████████████████████████████████] 100% [00m:00s] "D"'s SIDs validated
```

We also need to declare the exposures along the simulation.

``` r

# add SNV and indel default exposures. This will be used from simulated time 0
# up to the successive exposure change.
m_engine$add_exposure(coefficients = c(SBS13 = 0.2, SBS1 = 0.8))
m_engine$add_exposure(c(ID2 = 0.6, ID13 = 0.2, ID21 = 0.2))

# add a new SNV exposure that will be used from simulated
# time 100 up to the successive exposure change.
m_engine$add_exposure(time = 100, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))
m_engine$add_exposure(time = 120, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5,
                                    ID1 = 0.8, ID9 = 0.2))

m_engine
#> MutationEngine
#>  Passenger rates
#>    "A":
#>       [0,inf): {SNV: 2e-09, indel: 2e-09}
#>    "B":
#>       [0,inf): {SNV: 1e-09, indel: 2e-09, CNA: 1e-09}
#>    "C":
#>       [0,inf): {SNV: 2e-08, indel: 2e-09}
#>    "D":
#>       [0,inf): {SNV: 2e-09, indel: 2e-09}
#> 
#>  Driver mutations
#>    "A":
#>        (chr22(16085675)[GCCTCCCGA>G]) on random allele
#>        CNA("D",chr22(22453799), len: 200000, allele: 1)
#>    "B":
#>        (chr22(16050185)[A>C]) on allele 1
#>        CNA("A",chr22(16485130), len: 200000)
#>    "C":
#>        (chr22(32786322)[T>G]) on random allele
#>        DGCR8 P26L (chr22(20073563)[C>T]) on allele 1
#>    "D":
#>        (chr22(51240420)[G>T]) on random allele
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#>      [0, 100[: {"SBS1": 0.8, "SBS13": 0.2}
#>      [100, 120[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}
#>      [120, ∞[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}
#> 
#>    indel Timed Exposures
#>      [0, 120[: {"ID13": 0.2, "ID2": 0.6, "ID21": 0.2}
#>      [120, ∞[: {"ID1": 0.8, "ID9": 0.2}
```

We are now ready to build the phylogenetic forest.

``` r

# place mutations on the sample forest assuming 1000 pre-neoplastic SNVs and
# 500 indels
phylo_forest <- m_engine$place_mutations(sample_forest, 1000, 500)
#>  [█---------------------------------------] 0% [00m:00s] Placing mutations [████████████████████████████████████████] 100% [00m:00s] Mutations placed
```

The phylogenetic forest maintains the genome mutations (SBSs, indels,
and CNAs) of all the sampled cells.

``` r

mutations <- phylo_forest$get_sampled_cell_mutations()

head(mutations)
#>   chr     from allele ref alt cause         nature cell_id
#> 1  22 16066623      0   G   T  SBS1 pre-neoplastic 1009865
#> 2  22 16154367      0   G   C  SBS1 pre-neoplastic 1009865
#> 3  22 16230375      0   A   C  SBS1 pre-neoplastic 1009865
#> 4  22 16337416      0   T   C  SBS1 pre-neoplastic 1009865
#> 5  22 16460766      0   G   T  SBS1 pre-neoplastic 1009865
#> 6  22 16517421      0   C   G  SBS1 pre-neoplastic 1009865

CNAs <- phylo_forest$get_sampled_cell_CNAs()

head(CNAs)
#>   chr    begin      end type allele src.allele cause    nature cell_id
#> 1  22 22453799 22653798    D      1         NA          driver 1009865
#> 2  22 22453799 22653798    D      1         NA          driver 1004918
#> 3  22 16485130 16685129    A      2          1          driver 1004918
#> 4  22 50327093 50947147    A      3          1       passenger 1004918
#> 5  22 44105633 46776242    A      4          0       passenger 1004918
#> 6  22 22708735 24347001    A      5          0       passenger 1004918
```

## Sequencing level

We can simulate the sequencing of the collected samples by using the
phylogenetic forest. For each SBS and indel mutation in the phylogenetic
forest, the output reports the number of occurrences in the simulated
reads, the coverage of each mutation locus, and the corresponding VAF in
each of the samples.

``` r

# simulate the sequencing without normal sample and avoid progress bar
seq_results <- simulate_seq(phylo_forest, coverage = 30,
                            with_normal_sample = FALSE, quiet = TRUE)

head(seq_results$mutations)
#>   chr     from ref alt cause         nature S1.NV S1.DP     S1.VAF S2.NV S2.DP
#> 1  22 16050185   A   C     B         driver     0    29 0.00000000    15    29
#> 2  22 16053663   T   C  SBS2      passenger     1    27 0.03703704     0    29
#> 3  22 16057599   A   G  SBS2      passenger     1    25 0.04000000     0    38
#> 4  22 16060299   A   G  SBS3      passenger     1    29 0.03448276     0    28
#> 5  22 16066623   G   T  SBS1 pre-neoplastic    16    30 0.53333333    17    34
#> 6  22 16077397   A   C  SBS1 pre-neoplastic    12    23 0.52173913    26    39
#>      S2.VAF S3.NV S3.DP    S3.VAF
#> 1 0.5172414    11    21 0.5238095
#> 2 0.0000000     0    20 0.0000000
#> 3 0.0000000     0    14 0.0000000
#> 4 0.0000000     0    18 0.0000000
#> 5 0.5000000     9    19 0.4736842
#> 6 0.6666667    14    21 0.6666667
```

ProCESS allows to visualise the sequencing results.

``` r

plot_VAF_histogram(seq_results, labels = seq_results$mutations["nature"],
                   cuts = c(0.02, 1))
```

![The VAF histogram labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-19-1.png)

Marginal distributions can also be plotted.

``` r

plot_VAF_marginals(seq_results, labels = seq_results$mutations["nature"],
                   samples = c("S1", "S2", "S3"))
#> [[1]]
```

![The VAF marginal distribution S1 vs S2 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-20-1.png)

    #> 
    #> [[2]]

![The VAF marginal distribution S1 vs S3 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-20-2.png)

    #> 
    #> [[3]]

![The VAF marginal distribution S2 vs S3 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-20-3.png)

In S1 vs S3 figure, we can some spot passenger mutations occurring in
both S1 and S3. Let us identify these mutations.

``` r

seq_results$mutations %>% filter(nature == "passenger" &
                                   S1.VAF > 0 & S3.VAF > 0)
#>   chr     from ref alt cause    nature S1.NV S1.DP    S1.VAF S2.NV S2.DP S2.VAF
#> 1  22 26625316   T   C  SBS1 passenger    21    40 0.5250000     0    29      0
#> 2  22 33438401   A   G  SBS1 passenger    15    32 0.4687500     0    23      0
#> 3  22 33816566   A   C  SBS1 passenger    18    37 0.4864865     0    29      0
#> 4  22 42591172   T   C SBS13 passenger    15    29 0.5172414     0    27      0
#>   S3.NV S3.DP     S3.VAF
#> 1     1    66 0.01515152
#> 2     1    29 0.03448276
#> 3     2    26 0.07692308
#> 4     1    53 0.01886792
```
