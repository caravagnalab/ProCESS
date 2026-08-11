# A four-mutants example

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

We build a tumour with a single mutant `A` and grow it up to 10 cells.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` `` ``# set the seed of the random number generator for repeatability`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``0``)`` `` ``sim`` ``<-`` `[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)`(``width ``=`` ``600``, height ``=`` ``600``)`` `` ``# set the delta time between two time series samples`` ``sim``$``history_delta`` ``<-`` ``10`` `` ``# avoid drift`` ``sim``$``death_activation_level`` ``<-`` ``50`` `` ``sim``$``add_mutant``(``"A"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.12``, death ``=`` ``0.05``)``)`` `` ``sim``$``place_cell``(``"A"``, ``300``, ``300``)`` ``sim``$``run_up_to_size``(``"A"``, ``10``)`

We add a new mutant `B` as a sub-clone of `A` and let the simulation
evolve until `B` reaches 30 cells.

`sim``$``add_mutant``(``"B"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.145``, death ``=`` ``0.06``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"B"``)`` `` ``sim``$``run_up_to_size``(``"B"``, ``30``)`` ``#> [█████-----------------------------------] 10% [00m:00s] Saving snapshot [███████---------------------------------] 16% [00m:00s] Saving snapshot [██████████------------------------------] 23% [00m:00s] Saving snapshot [██████████████--------------------------] 33% [00m:00s] Saving snapshot [█████████████████████████---------------] 60% [00m:00s] Saving snapshot [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`

We add a further mutant `C` as a sub-clones for `A` and let the tumour
grow again up to the point at which `C` consists of 25000 cells.

`sim``$``add_mutant``(``"C"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.15``, death ``=`` ``0.06``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"C"``)`` `` ``sim``$``run_up_to_size``(``"C"``, ``25000``)`` ``#> [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 1% [00m:00s] Saving snapshot [█---------------------------------------] 1% [00m:00s] Saving snapshot [█---------------------------------------] 1% [00m:00s] Saving snapshot [█---------------------------------------] 1% [00m:00s] Saving snapshot [█---------------------------------------] 2% [00m:00s] Saving snapshot [█---------------------------------------] 2% [00m:00s] Saving snapshot [█---------------------------------------] 2% [00m:00s] Saving snapshot [██--------------------------------------] 3% [00m:00s] Saving snapshot [██--------------------------------------] 3% [00m:00s] Saving snapshot [██--------------------------------------] 4% [00m:00s] Saving snapshot [██--------------------------------------] 4% [00m:00s] Saving snapshot [███-------------------------------------] 5% [00m:00s] Saving snapshot [███-------------------------------------] 5% [00m:00s] Saving snapshot [███-------------------------------------] 6% [00m:00s] Saving snapshot [███-------------------------------------] 6% [00m:00s] Saving snapshot [███-------------------------------------] 7% [00m:00s] Saving snapshot [███-------------------------------------] 7% [00m:00s] Saving snapshot [████------------------------------------] 8% [00m:00s] Saving snapshot [████------------------------------------] 9% [00m:00s] Saving snapshot [████------------------------------------] 9% [00m:00s] Saving snapshot [█████-----------------------------------] 10% [00m:00s] Saving snapshot [█████-----------------------------------] 11% [00m:00s] Saving snapshot [█████-----------------------------------] 12% [00m:00s] Saving snapshot [█████-----------------------------------] 12% [00m:00s] Saving snapshot [██████----------------------------------] 13% [00m:00s] Saving snapshot [██████----------------------------------] 14% [00m:00s] Saving snapshot [███████---------------------------------] 15% [00m:00s] Saving snapshot [███████---------------------------------] 16% [00m:00s] Saving snapshot [███████---------------------------------] 17% [00m:00s] Saving snapshot [████████--------------------------------] 18% [00m:00s] Saving snapshot [████████--------------------------------] 19% [00m:00s] Saving snapshot [█████████-------------------------------] 20% [00m:00s] Saving snapshot [█████████-------------------------------] 21% [00m:00s] Saving snapshot [█████████-------------------------------] 22% [00m:00s] Saving snapshot [██████████------------------------------] 23% [00m:00s] Saving snapshot [██████████------------------------------] 24% [00m:00s] Saving snapshot [███████████-----------------------------] 25% [00m:00s] Saving snapshot [███████████-----------------------------] 26% [00m:00s] Cells: 14894 [███████████-----------------------------] 27% [00m:01s] Saving snapshot [████████████----------------------------] 28% [00m:01s] Saving snapshot [████████████----------------------------] 29% [00m:01s] Saving snapshot [█████████████---------------------------] 31% [00m:01s] Saving snapshot [█████████████---------------------------] 32% [00m:01s] Saving snapshot [██████████████--------------------------] 33% [00m:01s] Saving snapshot [███████████████-------------------------] 35% [00m:01s] Saving snapshot [███████████████-------------------------] 36% [00m:01s] Saving snapshot [███████████████-------------------------] 37% [00m:01s] Saving snapshot [████████████████------------------------] 39% [00m:01s] Saving snapshot [█████████████████-----------------------] 41% [00m:01s] Saving snapshot [█████████████████-----------------------] 42% [00m:01s] Saving snapshot [█████████████████-----------------------] 42% [00m:01s] Cells: 23326 [██████████████████----------------------] 44% [00m:02s] Saving snapshot [███████████████████---------------------] 45% [00m:02s] Saving snapshot [███████████████████---------------------] 47% [00m:02s] Saving snapshot [████████████████████--------------------] 49% [00m:02s] Saving snapshot [█████████████████████-------------------] 50% [00m:02s] Saving snapshot [█████████████████████-------------------] 52% [00m:02s] Saving snapshot [██████████████████████------------------] 54% [00m:02s] Saving snapshot [███████████████████████-----------------] 55% [00m:02s] Saving snapshot [███████████████████████-----------------] 56% [00m:02s] Cells: 30251 [███████████████████████-----------------] 57% [00m:03s] Saving snapshot [████████████████████████----------------] 58% [00m:03s] Saving snapshot [█████████████████████████---------------] 61% [00m:03s] Saving snapshot [██████████████████████████--------------] 63% [00m:03s] Saving snapshot [██████████████████████████--------------] 64% [00m:03s] Saving snapshot [███████████████████████████-------------] 66% [00m:03s] Saving snapshot [███████████████████████████-------------] 67% [00m:03s] Cells: 36251 [████████████████████████████------------] 68% [00m:04s] Saving snapshot [█████████████████████████████-----------] 70% [00m:04s] Saving snapshot [██████████████████████████████----------] 73% [00m:04s] Saving snapshot [███████████████████████████████---------] 75% [00m:04s] Saving snapshot [███████████████████████████████---------] 77% [00m:04s] Saving snapshot [███████████████████████████████---------] 77% [00m:04s] Cells: 41414 [████████████████████████████████--------] 79% [00m:05s] Saving snapshot [█████████████████████████████████-------] 81% [00m:05s] Saving snapshot [██████████████████████████████████------] 83% [00m:05s] Saving snapshot [███████████████████████████████████-----] 85% [00m:05s] Saving snapshot [███████████████████████████████████-----] 86% [00m:05s] Cells: 45718 [███████████████████████████████████-----] 87% [00m:06s] Saving snapshot [█████████████████████████████████████---] 90% [00m:06s] Saving snapshot [█████████████████████████████████████---] 92% [00m:06s] Saving snapshot [██████████████████████████████████████--] 94% [00m:06s] Saving snapshot [███████████████████████████████████████-] 95% [00m:06s] Cells: 50189 [███████████████████████████████████████-] 96% [00m:07s] Saving snapshot [████████████████████████████████████████] 98% [00m:07s] Saving snapshot [████████████████████████████████████████] 100% [00m:07s] Saving snapshot`

Then, we plot the tissue and the simulation state.

[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``sim``)`

![The tissue of a tumour having three mutants: A, B, and C. Both B, and
C are sub-clones of A. The three mutants have different grow rates and
appears at different simulated
time.](four_mutants_files/figure-html/unnamed-chunk-6-1.png)

[`plot_state`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_state.md)`(``sim``)`

![The current ratio of cells per species over all tumour
cells.](four_mutants_files/figure-html/unnamed-chunk-7-1.png)

## Tissue sampling level

We can collect two samples: “S1” and “S2”.

`sim``$``sample_cells``(``"S1"``, `[`c`](https://rdrr.io/r/base/c.html)`(``145``, ``230``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``215``, ``300``)``)`` ``sim``$``sample_cells``(``"S2"``, `[`c`](https://rdrr.io/r/base/c.html)`(``350``, ``300``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``420``, ``370``)``)`

We can plot the tissue simulation after the sampling. We label the
sampled region to improve readability, but it is not mandatory.

[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``sim``)`` ``+`` `` ``ggplot2``::`[`annotate`](https://ggplot2.tidyverse.org/reference/annotate.html)`(``"text"``, x ``=`` ``(``145`` ``+`` ``215``)`` ``/`` ``2``, y ``=`` ``(``230`` ``+`` ``300``)`` ``/`` ``2``,`` `` label ``=`` ``"S1"``, parse ``=`` ``TRUE``)`` ``+`` `` ``ggplot2``::`[`annotate`](https://ggplot2.tidyverse.org/reference/annotate.html)`(``"text"``, x ``=`` ``(``350`` ``+`` ``420``)`` ``/`` ``2``, y ``=`` ``(``300`` ``+`` ``370``)`` ``/`` ``2``,`` `` label ``=`` ``"S2"``, parse ``=`` ``TRUE``)`

![The three-mutants tumour just after collecting two samples: S1 and S2.
The former mainly contains C cells; the latter A
cells.](four_mutants_files/figure-html/unnamed-chunk-9-1.png)

After the sampling, we can add a new mutant `D` as a sub-clone of `B`,
let the simulation continue to evolve until the sum of the cardinalities
of mutants `C` and `D` is 100000, and sample the tissue again.

`sim``$``add_mutant``(``name ``=`` ``"D"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.8``, death ``=`` ``0.01``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"B"``)``, ``"D"``)`` `` ``sim``$``run_until``(``sim``$``var``(``"C"``)`` ``+`` ``sim``$``var``(``"D"``)`` ``==`` ``1e5``)`` ``#> [█---------------------------------------] 0% [00m:00s] Saving snapshot [█---------------------------------------] 1% [00m:00s] Saving snapshot [██--------------------------------------] 3% [00m:00s] Saving snapshot [██--------------------------------------] 4% [00m:00s] Cells: 51512 [██--------------------------------------] 4% [00m:00s] Saving snapshot [███-------------------------------------] 6% [00m:01s] Saving snapshot [████------------------------------------] 8% [00m:01s] Saving snapshot [█████-----------------------------------] 11% [00m:01s] Saving snapshot [█████-----------------------------------] 11% [00m:01s] Cells: 59553 [██████----------------------------------] 13% [00m:02s] Saving snapshot [███████---------------------------------] 16% [00m:02s] Saving snapshot [████████--------------------------------] 19% [00m:02s] Cells: 68461 [████████--------------------------------] 19% [00m:02s] Saving snapshot [██████████------------------------------] 23% [00m:03s] Saving snapshot [███████████-----------------------------] 27% [00m:03s] Saving snapshot [████████████----------------------------] 28% [00m:03s] Cells: 77101 [█████████████---------------------------] 31% [00m:04s] Saving snapshot [███████████████-------------------------] 36% [00m:04s] Saving snapshot [████████████████------------------------] 38% [00m:04s] Cells: 85762 [█████████████████-----------------------] 41% [00m:05s] Saving snapshot [███████████████████---------------------] 47% [00m:05s] Saving snapshot [████████████████████--------------------] 48% [00m:05s] Cells: 94769 [██████████████████████------------------] 52% [00m:06s] Saving snapshot [████████████████████████----------------] 58% [00m:06s] Cells: 104426 [████████████████████████----------------] 58% [00m:06s] Saving snapshot [███████████████████████████-------------] 66% [00m:07s] Saving snapshot [████████████████████████████------------] 69% [00m:07s] Cells: 113441 [██████████████████████████████----------] 73% [00m:08s] Saving snapshot [█████████████████████████████████-------] 81% [00m:08s] Cells: 122736 [█████████████████████████████████-------] 81% [00m:08s] Saving snapshot [█████████████████████████████████████---] 90% [00m:09s] Saving snapshot [█████████████████████████████████████---] 91% [00m:09s] Cells: 131817 [████████████████████████████████████████] 98% [00m:10s] Saving snapshot [████████████████████████████████████████] 100% [00m:10s] Saving snapshot`` `` ``sim``$``sample_cells``(``"S3"``, `[`c`](https://rdrr.io/r/base/c.html)`(``100``, ``400``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``170``, ``470``)``)`` `` `[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``sim``)`` ``+`` `` ``ggplot2``::`[`annotate`](https://ggplot2.tidyverse.org/reference/annotate.html)`(``"text"``, x ``=`` ``(``100`` ``+`` ``170``)`` ``/`` ``2``, y ``=`` ``(``400`` ``+`` ``470``)`` ``/`` ``2``,`` `` label ``=`` ``"S3"``, parse ``=`` ``TRUE``)`

![After collecting the samples S1 and S2, a more aggressive mutant D
arises from B and a further sample S3, containing both D and C cells, is
collected.](four_mutants_files/figure-html/unnamed-chunk-10-1.png)

` `[`plot_state`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_state.md)`(``sim``)`

![The simulation state after collecting
S3.](four_mutants_files/figure-html/unnamed-chunk-10-2.png)

The time-series plot represents the species cardinalities along the
simulation.

[`plot_timeseries`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_timeseries.md)`(``sim``)`

![The time series plot represent the species cardinalities along the
simulation.](four_mutants_files/figure-html/unnamed-chunk-11-1.png)

The Muller plot, instead, shows the percentage of cells in each species
over the course of the simulation.

[`plot_muller`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_muller.md)`(``sim``)`

![The Muller plot overviews the percentage of cells in each species
along the simulation
evolution.](four_mutants_files/figure-html/unnamed-chunk-12-1.png)

We can build the sample forest and plot it.

`sample_forest`` ``<-`` ``sim``$``get_sample_forest``(``)`` `` `[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`` ``#> `` ``#> Attaching package: 'dplyr'`` ``#> The following objects are masked from 'package:stats':`` ``#> `` ``#> filter, lag`` ``#> The following objects are masked from 'package:base':`` ``#> `` ``#> intersect, setdiff, setequal, union`` `` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``sample_forest``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``sample_forest``)`

![The sample forest of the collected
samples.](four_mutants_files/figure-html/unnamed-chunk-13-1.png)

## Cell genetic characterisation level

We first build the mutation engine to label each cell in the sample
forest. The mutation engine pre-defined setups automatically download
all the needed files. All the pre-configured setups, but “demo” requires
a COSMIC account to download the signature files from the COSMIC site.

`# build the mutation engine according to the pre-defined setup "demo"`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code``=``"demo"``,`` `` COSMIC_account ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``email ``=`` ``"foo@bar.com"``,`` `` password ``=`` ``"#########"``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading context index [████████████████████████████████████████] 100% [00m:00s] Context index loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Loading RS index [██████████████--------------------------] 33% [00m:01s] Loading RS index [███████████████████████████-------------] 66% [00m:02s] Loading RS index [████████████████████████████████████████] 98% [00m:03s] Loading RS index [████████████████████████████████████████] 100% [00m:03s] RS index loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`

We must genomically characterise the tissue simulation mutants,
providing passenger mutation rates and the list of driver mutations.

`# genetically characterise the mutants`` ``m_engine``$``add_mutant``(``"A"``, passenger_rates ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``2e-9``, indel ``=`` ``2e-9``)``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md)`(``"22"``, ``16085675``, ``"GCCTCCCGA"``, ``"G"``)``,`` `` `[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md)`(``"D"``, ``"22"``, ``22453799``, ``200000``,`` `` allele ``=`` ``1``)``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs [█---------------------------------------] 0% [00m:00s] Found 22 [█---------------------------------------] 0% [00m:00s] Reading 22 [█████████████---------------------------] 32% [00m:01s] Reading 22 [███████████████████████████-------------] 66% [00m:02s] Reading 22 [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated`` `` ``m_engine``$``add_mutant``(``"B"``, passenger_rates ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``1e-9``, indel ``=`` ``2e-9``,`` `` CNA ``=`` ``1e-9``)``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``16050185``, ``"C"``, allele ``=`` ``1``)``,`` `` `[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md)`(``type ``=`` ``"A"``, chr ``=`` ``"22"``,`` `` from ``=`` ``16485130``, len ``=`` ``200000``)``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated`` `` ``m_engine``$``add_mutant``(``"C"``, passenger_rates ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``SNV ``=`` ``2e-8``, indel ``=`` ``2e-9``)``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``32786322``, ``"G"``)``,`` `` `[`list`](https://rdrr.io/r/base/list.html)`(``"DGCR8 P26L"``, allele ``=`` ``1``)``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "C" SIDs [████████████████████████████████████████] 100% [00m:00s] "C"'s SIDs validated`` `` ``m_engine``$``add_mutant``(``"D"``, passenger_rates ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``SNV ``=`` ``2e-9``, indel ``=`` ``2e-9``)``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``51240420``, ``"T"``)``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "D" SIDs [████████████████████████████████████████] 100% [00m:00s] "D"'s SIDs validated`

We also need to declare the exposures along the simulation.

`# add SNV and indel default exposures. This will be used from simulated time 0`` ``# up to the successive exposure change.`` ``m_engine``$``add_exposure``(``coefficients ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SBS13 ``=`` ``0.2``, SBS1 ``=`` ``0.8``)``)`` ``m_engine``$``add_exposure``(`[`c`](https://rdrr.io/r/base/c.html)`(``ID2 ``=`` ``0.6``, ID13 ``=`` ``0.2``, ID21 ``=`` ``0.2``)``)`` `` ``# add a new SNV exposure that will be used from simulated`` ``# time 100 up to the successive exposure change.`` ``m_engine``$``add_exposure``(``time ``=`` ``100``, `[`c`](https://rdrr.io/r/base/c.html)`(``SBS5 ``=`` ``0.3``, SBS2 ``=`` ``0.2``, SBS3 ``=`` ``0.5``)``)`` ``m_engine``$``add_exposure``(``time ``=`` ``120``, `[`c`](https://rdrr.io/r/base/c.html)`(``SBS5 ``=`` ``0.3``, SBS2 ``=`` ``0.2``, SBS3 ``=`` ``0.5``,`` `` ID1 ``=`` ``0.8``, ID9 ``=`` ``0.2``)``)`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> "A":`` ``#> [0,inf): {SNV: 2e-09, indel: 2e-09}`` ``#> "B":`` ``#> [0,inf): {SNV: 1e-09, indel: 2e-09, CNA: 1e-09}`` ``#> "C":`` ``#> [0,inf): {SNV: 2e-08, indel: 2e-09}`` ``#> "D":`` ``#> [0,inf): {SNV: 2e-09, indel: 2e-09}`` ``#> `` ``#> Driver mutations`` ``#> "A":`` ``#> (chr22(16085675)[GCCTCCCGA>G]) on random allele`` ``#> CNA("D",chr22(22453799), len: 200000, allele: 1)`` ``#> "B":`` ``#> (chr22(16050185)[A>C]) on allele 1`` ``#> CNA("A",chr22(16485130), len: 200000)`` ``#> "C":`` ``#> (chr22(32786322)[T>G]) on random allele`` ``#> DGCR8 P26L (chr22(20073563)[C>T]) on allele 1`` ``#> "D":`` ``#> (chr22(51240420)[G>T]) on random allele`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> [0, 100[: {"SBS1": 0.8, "SBS13": 0.2}`` ``#> [100, 120[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}`` ``#> [120, ∞[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}`` ``#> `` ``#> indel Timed Exposures`` ``#> [0, 120[: {"ID13": 0.2, "ID2": 0.6, "ID21": 0.2}`` ``#> [120, ∞[: {"ID1": 0.8, "ID9": 0.2}`

We are now ready to build the phylogenetic forest.

`# place mutations on the sample forest assuming 1000 pre-neoplastic SNVs and`` ``# 500 indels`` ``phylo_forest`` ``<-`` ``m_engine``$``place_mutations``(``sample_forest``, ``1000``, ``500``)`` ``#> [█---------------------------------------] 0% [00m:00s] Placing mutations [████████████████████████████████████████] 100% [00m:00s] Mutations placed`

The phylogenetic forest maintains the genome mutations (SBSs, indels,
and CNAs) of all the sampled cells.

`mutations`` ``<-`` ``phylo_forest``$``get_sampled_cell_mutations``(``)`` `` `[`head`](https://rdrr.io/r/utils/head.html)`(``mutations``)`` ``#> chr from allele ref alt cause nature cell_id`` ``#> 1 22 16066623 0 G T SBS1 pre-neoplastic 1009865`` ``#> 2 22 16154367 0 G C SBS1 pre-neoplastic 1009865`` ``#> 3 22 16230375 0 A C SBS1 pre-neoplastic 1009865`` ``#> 4 22 16337416 0 T C SBS1 pre-neoplastic 1009865`` ``#> 5 22 16460766 0 G T SBS1 pre-neoplastic 1009865`` ``#> 6 22 16517421 0 C G SBS1 pre-neoplastic 1009865`` `` ``CNAs`` ``<-`` ``phylo_forest``$``get_sampled_cell_CNAs``(``)`` `` `[`head`](https://rdrr.io/r/utils/head.html)`(``CNAs``)`` ``#> chr begin end type allele src.allele cause nature cell_id`` ``#> 1 22 22453799 22653798 D 1 NA driver 1009865`` ``#> 2 22 22453799 22653798 D 1 NA driver 1385787`` ``#> 3 22 16485130 16685129 A 2 1 driver 1385787`` ``#> 4 22 50327093 50947147 A 3 1 passenger 1385787`` ``#> 5 22 44105633 46776242 A 4 0 passenger 1385787`` ``#> 6 22 22708735 24347001 A 5 0 passenger 1385787`

## Sequencing level

We can simulate the sequencing of the collected samples by using the
phylogenetic forest. For each SBS and indel mutation in the phylogenetic
forest, the output reports the number of occurrences in the simulated
reads, the coverage of each mutation locus, and the corresponding VAF in
each of the samples.

`# simulate the sequencing without normal sample and avoid progress bar`` ``seq_results`` ``<-`` `[`simulate_seq`](https://caravagnalab.github.io/ProCESS/1.3/reference/simulate_seq.md)`(``phylo_forest``, coverage ``=`` ``30``,`` `` with_normal_sample ``=`` ``FALSE``, quiet ``=`` ``TRUE``)`` `` `[`head`](https://rdrr.io/r/utils/head.html)`(``seq_results``$``mutations``)`` ``#> chr from ref alt cause nature S1.NV S1.DP S1.VAF S2.NV S2.DP`` ``#> 1 22 16050185 A C B driver 0 29 0.00000000 13 31`` ``#> 2 22 16065249 A G SBS3 passenger 3 33 0.09090909 0 25`` ``#> 3 22 16066623 G T SBS1 pre-neoplastic 16 30 0.53333333 19 39`` ``#> 4 22 16077397 A C SBS1 pre-neoplastic 12 23 0.52173913 14 33`` ``#> 5 22 16078428 T G SBS1 pre-neoplastic 17 39 0.43589744 17 29`` ``#> 6 22 16081360 C G SBS1 pre-neoplastic 17 32 0.53125000 14 32`` ``#> S2.VAF S3.NV S3.DP S3.VAF`` ``#> 1 0.4193548 20 30 0.6666667`` ``#> 2 0.0000000 0 25 0.0000000`` ``#> 3 0.4871795 15 24 0.6250000`` ``#> 4 0.4242424 13 34 0.3823529`` ``#> 5 0.5862069 6 24 0.2500000`` ``#> 6 0.4375000 9 21 0.4285714`

ProCESS allows to visualise the sequencing results.

[`plot_VAF_histogram`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF_histogram.md)`(``seq_results``, labels ``=`` ``seq_results``$``mutations``[``"nature"``]``,`` `` cuts ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.02``, ``1``)``)`

![The VAF histogram labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-20-1.png)

Marginal distributions can also be plotted.

[`plot_VAF_marginals`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_VAF_marginals.md)`(``seq_results``, labels ``=`` ``seq_results``$``mutations``[``"nature"``]``,`` `` samples ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"S1"``, ``"S2"``, ``"S3"``)``)`` ``#> [[1]]`

![The VAF marginal distribution S1 vs S2 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-21-1.png)

    #> 
    #> [[2]]

![The VAF marginal distribution S1 vs S3 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-21-2.png)

    #> 
    #> [[3]]

![The VAF marginal distribution S2 vs S3 labelled by mutation
class.](four_mutants_files/figure-html/unnamed-chunk-21-3.png)

In S1 vs S3 figure, we can some spot passenger mutations occurring in
both S1 and S3. Let us identify these mutations.

`seq_results``$``mutations`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``nature`` ``==`` ``"passenger"`` ``&`` `` ``S1.VAF`` ``>`` ``0`` ``&`` ``S3.VAF`` ``>`` ``0``)`` ``#> [1] chr from ref alt cause nature S1.NV S1.DP S1.VAF S2.NV `` ``#> [11] S2.DP S2.VAF S3.NV S3.DP S3.VAF`` ``#> <0 rows> (or 0-length row.names)`
