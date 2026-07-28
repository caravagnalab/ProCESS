# Infinite Site Model

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

The mutation engine places mutations on the sampled cell genome
according to the [infinite sites
model](https://en.wikipedia.org/wiki/Infinite_sites_model) by default.
In particular, any new mutation is placed on a locus whose context is
mutation-free.

Let us consider the phylogenetic forest as built in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)
and verify whether it satisfies the infinite site model.

[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`` ``#> `` ``#> Attaching package: 'dplyr'`` ``#> The following objects are masked from 'package:stats':`` ``#> `` ``#> filter, lag`` ``#> The following objects are masked from 'package:base':`` ``#> `` ``#> intersect, setdiff, setequal, union`` `` ``# this function verifies whether any mutation arose in`` ``# two unrelated cells represented in a phylogenetic forest`` ``test_infinite_sites_model`` ``<-`` ``function``(``phylo_forest``)`` ``{`` `` ``# extract non-germinal SNVs that appear multiple times either`` `` ``# in the same cell or different cells`` `` ``snvs`` ``<-`` ``phylo_forest``$``get_sampled_cell_mutations``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``nature`` ``!=`` ``"germinal"``,`` `` `[`nchar`](https://rdrr.io/r/base/nchar.html)`(``.data``$``ref``)`` ``==`` ``1``,`` `` `[`nchar`](https://rdrr.io/r/base/nchar.html)`(``.data``$``alt``)`` ``==`` ``1``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``.data``$``chr``, ``.data``$``from``, ``.data``$``ref``, ``.data``$``alt``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``n`` ``>`` ``1``)`` `` `` ``# search for an SNV that independently occurred in two unrelated cells`` `` ``first_occurrences`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``)`` `` ``row`` ``<-`` ``1`` `` ``while`` ``(`[`length`](https://rdrr.io/r/base/length.html)`(``first_occurrences``)`` ``<`` ``2`` ``&&`` ``row`` ``<=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``snvs``)``)`` ``{`` `` ``snv`` ``<-`` `[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``snvs``[``row``, ``"chr"``]``, ``snvs``[``row``, ``"from"``]``,`` `` ref ``=`` ``snvs``[``row``, ``"ref"``]``, alt ``=`` ``snvs``[``row``, ``"alt"``]``)`` `` `` ``first_occurrences`` ``<-`` ``phylo_forest``$``get_first_occurrences``(``snv``)`` `` ``row`` ``<-`` ``row`` ``+`` ``1`` `` ``}`` `` `` ``# if the last handled SNV independently occurred in two unrelated`` `` ``# cells at least`` `` ``if`` ``(`[`length`](https://rdrr.io/r/base/length.html)`(``first_occurrences``)`` ``>=`` ``2``)`` ``{`` `` `` ``# print a message containing the two cells`` `` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"SNV('"``, ``snv``$``get_chromosome``(``)``, ``"'',"``,`` `` ``snv``$``get_position_in_chromosome``(``)``, ``",'"``, ``snv``$``get_ref``(``)``,`` `` ``"','"``, ``snv``$``get_alt``(``)``,`` `` ``"') independently arises in cells "``, ``first_occurrences``[``1``]``,`` `` ``" and "``, ``first_occurrences``[``2``]``)`` `` ``}`` ``else`` ``{`` `` `[`print`](https://rdrr.io/r/base/print.html)`(``"Every mutation arises exclusively in one cell"``)`` `` ``}`` ``}`` `` ``# test whether the infinite sites conditions hold in the built forest`` ``test_infinite_sites_model``(``phylo_forest``)`` ``#> [1] "Every mutation arises exclusively in one cell"`

This behaviour can be changed by using the mutation engine property
[`MutationEngine$infinite_sites_model`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-infinite_sites_model.md).
This property is a Boolean flag that enables/disables the infinite sites
model.

`# establish whether the infinite sites model is used`` ``m_engine``$``infinite_sites_model`` ``#> [1] TRUE`` `` ``# disable it`` ``m_engine``$``infinite_sites_model`` ``<-`` ``FALSE`

When the infinite sites model is disabled,
[`MutationEngine$place_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-place_mutations.md)
may place two mutations in the same locus of different alleles of the
same genome or the same mutation in the same locus of one allele of the
genomes of two cells that are not each other ancestors.

`# test whether the infinite sites model is enable`` ``m_engine``$``infinite_sites_model`` ``#> [1] FALSE`` `` ``# place the mutations on the same sample forest used above`` ``phylo_forest2`` ``<-`` ``m_engine``$``place_mutations``(``sample_forest``, ``1000``, ``500``)`` ``#> [█---------------------------------------] 0% [00m:00s] Placing mutations [████████████████████████████████████████] 100% [00m:00s] Mutations placed`` `` ``# test whether the infinite sites conditions hold in the new forest`` ``test_infinite_sites_model``(``phylo_forest2``)`` ``#> [1] "SNV('22'',16144992,'C','T') independently arises in cells 114456 and 169972"`
