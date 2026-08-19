# Growth Models

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS/CLONES simulate the tumour growth according to two alternative
models: the “border” growth model and the homogeneous-growth model.

The former exclusively admits duplication of cells that have space to
duplicate (i.e., cells on the cancer border or near to dead cells). The
latter allows cell duplications everywhere in the tumour.

The “border” growth model is the default one. However, users can switch
to the homogeneous growth model by assigning the `TissueSimulation`’s
`border_growth_model` Boolean field. By setting it to `FALSE`, the
homogeneous-growth model is used. If, instead, it is set to `TRUE`, the
simulation evolves according to the “border” growth model.

In the rest of this article, we clarify the differences among the
supported growth models by showing that they produce distinct evolutions
of the same cancer model.

### Homogeneous Growth

First of all, we need to create a new `TissueSimulation` object and set
`border_growth_model` to `FALSE`. The simulation has two epigenetic
states: E1 and E2.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` ``#> `` ``#> Attaching package: 'ProCESS'`` ``#> The following object is masked from 'package:utils':`` ``#> `` ``#> example`` `` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``0``)`` `` ``sim`` ``<-`` `[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)`(``"Homogeneous Growth"``, epigenetic_states ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"E1"``, ``"E2"``)``)`` `` ``# set the homogeneous growth model`` ``sim``$``border_growth_model`` ``<-`` ``FALSE`` `` ``# set the death activation level to avoid drift`` ``sim``$``death_activation_level`` ``<-`` ``50`

Add a mutant `A`, set its simulation rates, and let the simulation
evolve until there are 1300 cells of species `A[E1]`. Then, collect two
samples, and let the simulation evolve again for 15 time units.

`# add a mutant`` ``sim``$``add_mutant``(``"A"``, `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.1``, death ``=`` ``0.08``,`` `` E2 ``=`` ``0.01``)``,`` `` ``"E2"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.1``, death ``=`` ``0.01``,`` `` E1 ``=`` ``0.01``)``)``)`` `` ``sim``$``place_cell``(``"A[E1]"``, ``500``, ``500``)`` `` ``# let the simulation evolve until "A[E1]" consists of 1300 cells`` ``sim``$``run_up_to_size``(``"A[E1]"``, ``1300``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``bbox_width`` ``<-`` ``15`` `` ``# takes two samples`` ``sim``$``sample_cells``(``"S_1_1"``,`` `` bottom_left ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``480``, ``480``)``,`` `` top_right ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``480`` ``+`` ``bbox_width``, ``480`` ``+`` ``bbox_width``)``)`` `` ``sim``$``sample_cells``(``"S_1_2"``,`` `` bottom_left ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``500``, ``500``)``,`` `` top_right ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``500`` ``+`` ``bbox_width``, ``500`` ``+`` ``bbox_width``)``)`` `` ``# let the simulation evolve again for 15 time units`` ``sim``$``run_up_to_time``(``sim``$``get_clock``(``)`` ``+`` ``15``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`

Add a new mutant `B`, let one of the cells in `A` generate a cell in
`B`, let the simulation evolve until there are 5000 cells in `B[E1]`.
Then, take two more samples.

`sim``$``add_mutant``(``"B"``, `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.8``, death ``=`` ``0.05``,`` `` E2 ``=`` ``0.05``)``,`` `` ``"E2"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.3``, death ``=`` ``0.05``,`` `` E1 ``=`` ``0.1``)``)``)`` `` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"B"``)`` `` ``sim``$``run_up_to_size``(``"B[E1]"``, ``5000``)`` ``#> [█---------------------------------------] 1% [00m:00s] Cells: 124738 [█████████████---------------------------] 31% [00m:01s] Cells: 179091 [████████████████████████████████████████] 100% [00m:01s] Saving snapshot`` `` ``ncells`` ``<-`` ``0.9`` ``*`` ``bbox_width`` ``*`` ``bbox_width`` `` ``bbox`` ``<-`` ``sim``$``search_sample``(`[`c`](https://rdrr.io/r/base/c.html)`(``"B"`` ``=`` ``ncells``)``, ``bbox_width``, ``bbox_width``)`` ``sim``$``sample_cells``(``"S_2_1"``, ``bbox``$``lower_corner``, ``bbox``$``upper_corner``)`` `` ``bbox`` ``<-`` ``sim``$``search_sample``(`[`c`](https://rdrr.io/r/base/c.html)`(``"A"`` ``=`` ``ncells``)``, ``bbox_width``, ``bbox_width``)`` ``sim``$``sample_cells``(``"S_2_2"``, ``bbox``$``lower_corner``, ``bbox``$``upper_corner``)`

Let us plot the tissue and the simulation Muller plot.

[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``sim``, num_of_bins ``=`` ``500``)`

![The tissue plot of a tumour grown according the uniform growth
model.](growth_models_files/figure-html/unnamed-chunk-5-1.png)

[`plot_muller`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_muller.md)`(``sim``)`

![The Muller plot of a tumour grown according the uniform growth
model.](growth_models_files/figure-html/unnamed-chunk-5-2.png)

Finally, let us build the sample forest.

[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`` ``#> `` ``#> Attaching package: 'dplyr'`` ``#> The following objects are masked from 'package:stats':`` ``#> `` ``#> filter, lag`` ``#> The following objects are masked from 'package:base':`` ``#> `` ``#> intersect, setdiff, setequal, union`` `` ``forest`` ``<-`` ``sim``$``get_sample_forest``(``)`` `` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of a tumour grown according the uniform growth
model.](growth_models_files/figure-html/unnamed-chunk-6-1.png)

We use a special `highlight` parameter to add the edges connecting cells
in each sample.

[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_1_1"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of a tumour grown according the uniform growth
model. The sub-forest of sample S_1_1 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-7-1.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_1_2"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of a tumour grown according the uniform growth
model. The sub-forest of sample S_1_2 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-7-2.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_2_1"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of a tumour grown according the uniform growth
model. The sub-forest of sample S_2_1 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-7-3.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_2_2"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of a tumour grown according the uniform growth
model. The sub-forest of sample S_2_2 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-7-4.png)

### “Border” Growth

New let us run the very same simulation by using the “border” growth
model.

[`set.seed`](https://rdrr.io/r/base/Random.html)`(``0``)`` `` ``sim`` ``<-`` `[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)`(``"Border Growth"``, epigenetic_states ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"E1"``, ``"E2"``)``)`` `` ``# setting the "border" growth model is not needed as the`` ``# border growth model is the default.`` ``sim``$``border_growth_model`` ``#> [1] TRUE`` `` ``# set the death activation level to avoid drift`` ``sim``$``death_activation_level`` ``<-`` ``50`

Let us repeat what we did in the homogeneous growth model example.

`# add a mutant`` ``sim``$``add_mutant``(``"A"``, `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.1``, death ``=`` ``0.08``,`` `` E2 ``=`` ``0.01``)``,`` `` ``"E2"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.1``, death ``=`` ``0.01``,`` `` E1 ``=`` ``0.01``)``)``)`` `` ``sim``$``place_cell``(``"A[E1]"``, ``500``, ``500``)`` `` ``# let the simulation evolve until "A[E1]" consists of 1300 cells`` ``sim``$``run_up_to_size``(``"A[E1]"``, ``1300``)`` ``#> [██████████████████████████--------------] 63% [00m:00s] Cells: 26181 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``bbox_width`` ``<-`` ``15`` `` ``# takes two samples`` ``sim``$``sample_cells``(``"S_1_1"``,`` `` bottom_left ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``480``, ``480``)``,`` `` top_right ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``480`` ``+`` ``bbox_width``, ``480`` ``+`` ``bbox_width``)``)`` `` ``sim``$``sample_cells``(``"S_1_2"``,`` `` bottom_left ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``500``, ``500``)``,`` `` top_right ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``500`` ``+`` ``bbox_width``, ``500`` ``+`` ``bbox_width``)``)`` `` ``# let the simulation evolve again for 15 time units`` ``sim``$``run_up_to_time``(``sim``$``get_clock``(``)`` ``+`` ``15``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``# add a new mutant`` ``sim``$``add_mutant``(``"B"``, `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.8``, death ``=`` ``0.05``,`` `` E2 ``=`` ``0.05``)``,`` `` ``"E2"`` ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``duplication ``=`` ``0.3``, death ``=`` ``0.05``,`` `` E1 ``=`` ``0.1``)``)``)`` `` ``# let one of the "A" cells generate a cell in "B"`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"B"``)`` `` ``# let the simulation evolve until "B[E1]" consists of 5000 cells`` ``sim``$``run_up_to_size``(``"B[E1]"``, ``5000``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``ncells`` ``<-`` ``0.9`` ``*`` ``bbox_width`` ``*`` ``bbox_width`` `` ``bbox`` ``<-`` ``sim``$``search_sample``(`[`c`](https://rdrr.io/r/base/c.html)`(``"B"`` ``=`` ``ncells``)``, ``bbox_width``, ``bbox_width``)`` ``sim``$``sample_cells``(``"S_2_1"``, ``bbox``$``lower_corner``, ``bbox``$``upper_corner``)`` `` ``bbox`` ``<-`` ``sim``$``search_sample``(`[`c`](https://rdrr.io/r/base/c.html)`(``"A"`` ``=`` ``ncells``)``, ``bbox_width``, ``bbox_width``)`` ``sim``$``sample_cells``(``"S_2_2"``, ``bbox``$``lower_corner``, ``bbox``$``upper_corner``)`

Let us have a look at the simulated tissue and plot the simulation
Muller plot.

[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``sim``, num_of_bins ``=`` ``500``)`

![The tissue plot of the above tumour grown according the border-driven
growth model.](growth_models_files/figure-html/unnamed-chunk-10-1.png)

[`plot_muller`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_muller.md)`(``sim``)`

![The Muller plot of the above tumour grown according the border-driven
growth model.](growth_models_files/figure-html/unnamed-chunk-10-2.png)

Once more, let us build the ancestor forest of the samples.

`forest`` ``<-`` ``sim``$``get_sample_forest``(``)`` `` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of the tumour grown according the border-driven
growth model.](growth_models_files/figure-html/unnamed-chunk-11-1.png)

We use a special `highlight` parameter to add the edges connecting cells
in each sample.

[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_1_1"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of the tumour grown according the border-driven
growth model. The sub-forest of sample S_1_1 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-12-1.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_1_2"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of the tumour grown according the border-driven
growth model. The sub-forest of sample S_1_2 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-12-2.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_2_1"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of the tumour grown according the border-driven
growth model. The sub-forest of sample S_2_1 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-12-3.png)

` `[`plot_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_forest.md)`(``forest``, highlight ``=`` ``"S_2_2"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`annotate_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/annotate_forest.md)`(``forest``)`

![The sample forest of the tumour grown according the border-driven
growth model. The sub-forest of sample S_2_2 is
highlighted.](growth_models_files/figure-html/unnamed-chunk-12-4.png)

It is easy to see the differences in the analogous plots of the above
examples.

First of all, the cancer growth is slower when subject to the “border”
growth model than when the homogeneous-growth model is used. This is
because the internal cells cannot duplicate in the former growth model;
thus, the number of cells active for duplication is always greater in it
than in the homogeneous-growth model.

Moreover, spatial closeness and closeness in the sample forest are
strictly related in the “border” growth model, whereas they appear to be
loosely related in the homogeneous growth model. These features can be
easily spotted in the forest plots when the sample `S_2_1` is selected
near the external tumour border.
