# Partitioning samples (FACS)

> *Note:* This article presents advanced topics on sequencing
> simulation. Refer to [this
> article](https://caravagnalab.github.io/ProCESS/1.3/articles/sequencing.md)
> for an introduction on the subject.

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS allows to partition samples according to cell features as
[FACS](https://www.sinobiological.com/category/fcm-facs-facs) does. The
method
[`PhylogeneticForest$partition_samples()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-partition_samples.md)
takes a labelling function, which accepts as input an object of the
class `SampledCell` and returns a string, and applies it to the sampled
cells. Then, it partitions the forest samples by using the obtained
labels. The resulting sub-samples are named according the format
`<sample name>_<cell label>`.

The method
[`PhylogeneticForest$partition_samples()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-partition_samples.md)
permanently changes the phylogenetic forest from which is called. If
users want to preserve the original forest, they must take care of it,
for instance by copying or saving it in a file.

Let us assume to have built the phylogenetic forest `phylo_forest` as
detailed in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md).
Information about the samples from which the forest has been built can
be obtained by using the
[`PhylogeneticForest$get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_samples_info.md).

``` r

phylo_forest$get_samples_info()
#>    name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 S_1_1  0  480  480  530  530         2576                 2576 445.2994
#> 2 S_1_2  1  500  500  550  550         1609                 1609 445.2994
#> 3 S_2_1  2  399  540  423  564          620                  620 506.0268
#> 4 S_2_2  3  549  365  573  389          620                  620 506.0268
#>   DNA_quantity equivalent_normal_cells
#> 1 528968984885                5155.184
#> 2 330460533638                3220.576
#> 3 128456775117                1251.904
#> 4 127235323680                1240.000
```

The forest samples can be split by epigenetic state as it follows.

``` r

# the labelling function parameter has type `SampledCell`
epi_labelling <- function(cell) {
  cell$epistate
}

phylo_forest$partition_samples(epi_labelling)

phylo_forest$get_samples_info()
#>       name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 S_1_1_E1  8  480  480  530  530           44                 2576 445.2994
#> 2 S_1_1_E2  4  480  480  530  530         2532                 2576 445.2994
#> 3 S_1_2_E1  7  500  500  550  550           28                 1609 445.2994
#> 4 S_1_2_E2  5  500  500  550  550         1581                 1609 445.2994
#> 5 S_2_1_E1  9  399  540  423  564           12                  620 506.0268
#> 6 S_2_1_E2  6  399  540  423  564          608                  620 506.0268
#> 7 S_2_2_E1 11  549  365  573  389          547                  620 506.0268
#> 8 S_2_2_E2 10  549  365  573  389           73                  620 506.0268
#>   DNA_quantity equivalent_normal_cells
#> 1   9031332253                88.01685
#> 2 519937652632              5067.16744
#> 3   5758351651                56.11929
#> 4 324702181987              3164.45696
#> 5   2493581895                24.30175
#> 6 125963193222              1227.60217
#> 7 112254390408              1094.00000
#> 8  14980933272               146.00000
```

The same method can be easily applied to partition sampled cells by
mutant name…

``` r

# loading again the original forest
phylo_forest <- load_phylogenetic_forest("phylo_forest.sff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# a mutant-name-based labelling function
mutant_labelling <- function(cell) {
  cell$mutant
}

phylo_forest$partition_samples(mutant_labelling)

phylo_forest$get_samples_info()
#>      name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox     time
#> 1 S_1_1_A 16  480  480  530  530         2576                 2576 445.2994
#> 2 S_1_2_A 17  500  500  550  550         1609                 1609 445.2994
#> 3 S_2_1_A 18  399  540  423  564          620                  620 506.0268
#> 4 S_2_2_A 19  549  365  573  389           11                  620 506.0268
#> 5 S_2_2_B 20  549  365  573  389          609                  620 506.0268
#>   DNA_quantity equivalent_normal_cells
#> 1 528968984885                5155.184
#> 2 330460533638                3220.576
#> 3 128456775117                1251.904
#> 4   2257400904                  22.000
#> 5 124977922776                1218.000
```

…, birth time…

``` r

# loading again the original forest
phylo_forest <- load_phylogenetic_forest("phylo_forest.sff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# a birth-time-based labelling function
birth_time_labelling <- function(cell) {
  if (cell$birth_time > 421) {
    return("YOUNG")
  }

  if (cell$birth_time > 321) {
    return("MIDDLE_AGED")
  }

  return("OLD")
}

phylo_forest$partition_samples(birth_time_labelling)

phylo_forest$get_samples_info()
#>                name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox
#> 1 S_1_1_MIDDLE_AGED 27  480  480  530  530         1342                 2576
#> 2         S_1_1_OLD 25  480  480  530  530          178                 2576
#> 3       S_1_1_YOUNG 30  480  480  530  530         1056                 2576
#> 4 S_1_2_MIDDLE_AGED 28  500  500  550  550          856                 1609
#> 5         S_1_2_OLD 26  500  500  550  550          106                 1609
#> 6       S_1_2_YOUNG 32  500  500  550  550          647                 1609
#> 7 S_2_1_MIDDLE_AGED 29  399  540  423  564           22                  620
#> 8       S_2_1_YOUNG 31  399  540  423  564          598                  620
#> 9       S_2_2_YOUNG 33  549  365  573  389          620                  620
#>       time DNA_quantity equivalent_normal_cells
#> 1 445.2994 275352014303              2683.50398
#> 2 445.2994  36611777738               356.80818
#> 3 445.2994 217005192844              2114.87212
#> 4 445.2994 175840546579              1713.69295
#> 5 445.2994  21761051636               212.07714
#> 6 445.2994 132858935423              1294.80615
#> 7 506.0268   4532003323                44.16764
#> 8 506.0268 123924771794              1207.73628
#> 9 506.0268 127235323680              1240.00000
```

…, mutations…

``` r

# loading again the original forest
phylo_forest <- load_phylogenetic_forest("phylo_forest.sff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# collect all the sample passenger indels
passenger_indels <- phylo_forest$get_sampled_cell_mutations() %>%
  filter(class == "passenger", type == "indel")

# get one of the passenger indels
p_indel <- passenger_indels[sample(seq_len(nrow(passenger_indels)), 1), ]

p_indel
#>      cell_id chr     from allele  ref alt  type cause     class
#> 5293  173106  22 35124166      0 AAGG   A indel  ID13 passenger

# a mutation-based labelling function that discriminates sampled cells
# containing indel `p_indel` from the other sampled cells
mutations_labelling <- function(cell) {
  has_indel <- nrow(cell$mutations %>%
                      filter(chr == p_indel[["chr"]],
                             from == p_indel[["from"]],
                             ref == p_indel[["ref"]],
                             alt == p_indel[["alt"]])) > 0

  if (has_indel) {
    return("HAS_MUTATION");
  }

  return("MISSES_MUTATION")
}

phylo_forest$partition_samples(mutations_labelling)

phylo_forest$get_samples_info()
#>                    name id xmin ymin xmax ymax tumour_cells
#> 1    S_1_1_HAS_MUTATION 38  480  480  530  530          660
#> 2 S_1_1_MISSES_MUTATION 39  480  480  530  530         1916
#> 3    S_1_2_HAS_MUTATION 40  500  500  550  550         1284
#> 4 S_1_2_MISSES_MUTATION 41  500  500  550  550          325
#> 5 S_2_1_MISSES_MUTATION 42  399  540  423  564          620
#> 6 S_2_2_MISSES_MUTATION 43  549  365  573  389          620
#>   tumour_cells_in_bbox     time DNA_quantity equivalent_normal_cells
#> 1                 2576 445.2994 135485896128               1320.4078
#> 2                 2576 445.2994 393483088757               3834.7765
#> 3                 1609 445.2994 263763522328               2570.5658
#> 4                 1609 445.2994  66697011310                650.0105
#> 5                  620 506.0268 128456775117               1251.9039
#> 6                  620 506.0268 127235323680               1240.0000
```

… or combination of these properties.
