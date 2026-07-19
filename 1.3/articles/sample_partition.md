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
class
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
and returns a string, and applies it to the sampled cells. Then, it
partitions the forest samples by using the obtained labels. The
resulting sub-samples are named according the format
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
#> 1 531196062192                5176.889
#> 2 330425908811                3220.239
#> 3 127222173306                1239.872
#> 4 127235323680                1240.000
```

The forest samples can be split by epigenetic state as it follows.

``` r

# the labelling function parameter has type `PhylogeneticForestNode`
epi_labelling <- function(forest_node) {
  forest_node$epistate_name
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
#> 1   9074081066                88.43346
#> 2 522121981126              5088.45530
#> 3   5751057816                56.04821
#> 4 324674850995              3164.19060
#> 5   2462619168                24.00000
#> 6 124759554138              1215.87184
#> 7 112254390408              1094.00000
#> 8  14980933272               146.00000
```

The same method can be easily applied to partition sampled cells by
mutant name…

``` r

# loading again the original forest
phylo_forest <- load_forest("phylo_forest.pff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# a mutant-name-based labelling function
mutant_labelling <- function(forest_node) {
  forest_node$mutant_name
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
#> 1 531196062192                5176.889
#> 2 330425908811                3220.239
#> 3 127222173306                1239.872
#> 4   2257400904                  22.000
#> 5 124977922776                1218.000
```

…, birth time…

``` r

# loading again the original forest
phylo_forest <- load_forest("phylo_forest.pff")
#>  [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded

# a birth-time-based labelling function
birth_time_labelling <- function(forest_node) {
  if (forest_node$birth_time > 421) {
    return("YOUNG")
  }

  if (forest_node$birth_time > 321) {
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
#> 1 445.2994 276897757834               2698.5684
#> 2 445.2994  36675607392                357.4302
#> 3 445.2994 217622696966               2120.8901
#> 4 445.2994 175779746283               1713.1004
#> 5 445.2994  21764677640                212.1125
#> 6 445.2994 132881484888               1295.0259
#> 7 506.0268   4514801808                 44.0000
#> 8 506.0268 122707371498               1195.8718
#> 9 506.0268 127235323680               1240.0000
```

…, mutations…

``` r

# loading again the original forest
phylo_forest <- load_forest("phylo_forest.pff")
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
  filter(nature == "passenger", nchar(ref) != nchar(alt))

# get one of the passenger indels and build the mutation
p_indel <- passenger_indels[sample(seq_along(passenger_indels), 1), ]

mutation <- Mutation(chr = p_indel[["chr"]], from = p_indel[["from"]],
                     ref = p_indel[["ref"]], alt = p_indel[["alt"]])

mutation
#> indel(chr: 22, from: 46070274, allele: random, ref: T, alt: TTG)

# find the nodes where the mutation occurs for the first time
first_node_ids <- phylo_forest$get_first_occurrences(mutation)

# a mutation-based labelling function that discriminates sampled cells
# that are descendant of the first nodes in which the mutation occurs, i.e.,
# sampled cells containing indel `p_indel`, from the other sampled cells
mutations_labelling <- function(forest_node) {

  while (forest_node$cell_id < first_node_ids) {
    forest_node <- forest_node$parent
  }

  if (forest_node$cell_id == first_node_ids) {
    return("HAS_MUTATION");
  }

  return("MISSES_MUTATION")
}

phylo_forest$partition_samples(mutations_labelling)

phylo_forest$get_samples_info()
#>                    name id xmin ymin xmax ymax tumour_cells
#> 1 S_1_1_MISSES_MUTATION 38  480  480  530  530         2576
#> 2 S_1_2_MISSES_MUTATION 39  500  500  550  550         1609
#> 3 S_2_1_MISSES_MUTATION 40  399  540  423  564          620
#> 4 S_2_2_MISSES_MUTATION 41  549  365  573  389          620
#>   tumour_cells_in_bbox     time DNA_quantity equivalent_normal_cells
#> 1                 2576 445.2994 531196062192                5176.889
#> 2                 1609 445.2994 330425908811                3220.239
#> 3                  620 506.0268 127222173306                1239.872
#> 4                  620 506.0268 127235323680                1240.000
```

… or combination of these properties.
