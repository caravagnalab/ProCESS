# Producing evolution videos

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

This article exclusively deals with producing videos of the tumour
evolution. Please, refer to [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/tissue_simulation.md)
for tumour evolution modelling.

## Collecting data

In order to produce a video, we need to collect simulation snapshots
during tumour evolution.
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)
performs this task by default at the end of
[`place_cell()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-place_cell.md),
[`sample_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-sample_cells.md),
[`mutate_progeny()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-mutate_progeny.md),
and `run_*()` calls.

Let us build a simulation and use
[`get_snapshot_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_snapshot_info.md)
to inspect the collected snapshots.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` ``#> `` ``#> Attaching package: 'ProCESS'`` ``#> The following object is masked from 'package:utils':`` ``#> `` ``#> example`` `` ``# set the seed of the random number generator for repeatability`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``0``)`` `` ``sim`` ``<-`` `[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)`(``width ``=`` ``600``, height ``=`` ``600``)`` `` ``# avoid drift`` ``sim``$``death_activation_level`` ``<-`` ``50`` `` ``# add a mutant`` ``sim``$``add_mutant``(``"A"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.12``, death ``=`` ``0.05``)``)`` `` ``# inspect the collected snapshots`` ``sim``$``get_snapshot_info``(``)`` ``#> [1] time clock cells file `` ``#> <0 rows> (or 0-length row.names)`` `` ``# place a cell and trigger the first snapshot`` ``sim``$``place_cell``(``"A"``, ``300``, ``300``)`` `` ``# count the collected snapshots`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` ``#> [1] 1`` `` ``# let the model evolve and trigger the second snapshot`` ``sim``$``run_up_to_size``(``"A"``, ``10``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``# count the collected snapshots`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` ``#> [1] 2`` `` ``# add a further mutant`` ``sim``$``add_mutant``(``"B"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.145``, death ``=`` ``0.06``)``)`` `` ``# one the "A" cell mutate into "B" and the third snapshot`` ``# is collected`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"B"``)`` `` ``# count the collected snapshots`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` ``#> [1] 3`` `` ``# let the model evolve and trigger the fourth snapshot`` ``sim``$``run_up_to_size``(``"B"``, ``30``)`` ``#> [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``# count the collected snapshots`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` ``#> [1] 4`` `` ``# collect a sample and get the fifth snapshot`` ``sim``$``sample_cells``(``"S1"``, `[`c`](https://rdrr.io/r/base/c.html)`(``145``, ``230``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``215``, ``300``)``)`` `` ``# count the collected snapshots`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` ``#> [1] 5`

The function
[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md),
which is exclusively available when the package
[av](https://CRAN.R-project.org/package=av) is installed, uses the
collected snapshot to build a video whose frames depicts the tissue as
it was when snapshots were collected.

The snippet

[`build_snapshot_video`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)`(``sim``, quiet ``=`` ``TRUE``)`

produces the following video.

## Producing more frames

Whenever the snapshots collected by default are not enough to completely
represent the tumour evolution, user can set custom snapshot triggers.
ProCESS supports three types of triggers: - computation time triggers,
e.g., collect a snapshot at least every 5 minutes - simulated time
triggers, e.g., collect a snapshot at least every 5 simulated time
units - number of tumour cells triggers, e.g., collect a snapshot at
least every time the overall number of tumour cells increased by 1000
cells.

The
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)
property
[`snapshot_triggers`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_triggers.md)
reports the current custom triggers and can be used to set them.

`# see the current custom snapshot triggers`` ``sim``$``snapshot_triggers`` ``#> list()`` `` ``# take a new snapshot every 1000 new tumour cells`` ``sim``$``snapshot_triggers`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``"number of cells"`` ``=`` ``1000``)`` `` ``# get the current simulated time`` ``sim``$``get_clock``(``)`` ``#> [1] 48.30854`` `` ``# count the snapshots collected before trigger changes`` ``previous_snapshots`` ``<-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` `` ``# let the model evolve until "B" consists of less than 2000 cells`` ``sim``$``run_up_to_size``(``"B"``, ``2000``)`` ``#> [██████████------------------------------] 23% [00m:00s] Saving snapshot [███████████████████---------------------] 46% [00m:00s] Saving snapshot [█████████████████████████████-----------] 70% [00m:00s] Saving snapshot [███████████████████████████████████████-] 95% [00m:00s] Saving snapshot [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``# get the simulated time at each snapshot`` ``snapshot_info`` ``<-`` ``sim``$``get_snapshot_info``(``)`` `` ``# show the number of tumour cells in each snapshot since the last default`` ``# trigger`` ``snapshot_info``$``cells``[``previous_snapshots``:`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``snapshot_info``)``]`` ``#> [1] 127 1127 2127 3127 4127 4304`

[`snapshot_triggers`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_triggers.md)
can also be used to reset custom triggers.

`# see the current custom snapshot triggers`` ``sim``$``snapshot_triggers`` ``` #> $`number of cells` ``` ``#> [1] 1000`` `` ``# take a new snapshot every 1000 new tumour cells`` ``sim``$``snapshot_triggers`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``"number of cells"`` ``=`` ``NULL``, ``"clock interval"`` ``=`` ``5``)`` `` ``# count the snapshots collected before trigger changes`` ``previous_snapshots`` ``<-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim``$``get_snapshot_info``(``)``)`` `` ``# let the model evolve for further 20 simulated time units`` ``sim``$``run_up_to_time``(``sim``$``get_clock``(``)`` ``+`` ``20``)`` ``#> [█████████████████████████████████████---] 92% [00m:00s] Saving snapshot [███████████████████████████████████████-] 95% [00m:00s] Saving snapshot [███████████████████████████████████████-] 97% [00m:00s] Saving snapshot [████████████████████████████████████████] 100% [00m:00s] Saving snapshot`` `` ``# get the simulated time at each snapshot`` ``snapshot_info`` ``<-`` ``sim``$``get_snapshot_info``(``)`` `` ``# show the simulated clock of each snapshot since trigger changes`` ``snapshot_info``$``clock``[``previous_snapshots``:`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``snapshot_info``)``]`` ``#> [1] 180.9084 185.9117 190.9132 195.9168 200.9110`

Custom snapshot triggers allow us to produce a more fluid evolution
video. For instance, let us consider the tumour model detailed in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/four_mutants.md)
and let us customise the snapshot triggers just after creating the
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)
object as follows.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` `` ``# set the seed of the random number generator for repeatability`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``0``)`` `` ``sim`` ``<-`` `[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)`(``width ``=`` ``600``, height ``=`` ``600``)`` `` ``# collect a snapshot every 5 simulated time units`` ``sim``$``snapshot_triggers`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``"clock interval"`` ``=`` ``5``)`

The remaining part of the model specification is the same detailed in
[this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/four_mutants.md).

`# avoid drift`` ``sim``$``death_activation_level`` ``<-`` ``50`` `` ``# add the mutant A`` ``sim``$``add_mutant``(``"A"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.12``, death ``=`` ``0.05``)``)`` `` ``# place a cell in the tissue and simulate it until 10 cells`` ``sim``$``place_cell``(``"A"``, ``300``, ``300``)`` ``sim``$``run_up_to_size``(``"A"``, ``10``, quiet ``=`` ``TRUE``)`` `` ``# add the mutant B and let mutate a border cell of A in B`` ``sim``$``add_mutant``(``"B"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.145``, death ``=`` ``0.06``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"B"``)`` `` ``# simulate the tissue up to 30 cells in B`` ``sim``$``run_up_to_size``(``"B"``, ``30``, quiet ``=`` ``TRUE``)`` `` ``# add the third mutant and let one cell of A mutate into C`` ``sim``$``add_mutant``(``"C"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.15``, death ``=`` ``0.06``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"A"``)``, ``"C"``)`` `` ``# simulate the tissue until C consists of 25000 cells`` ``sim``$``run_up_to_size``(``"C"``, ``25000``, quiet ``=`` ``TRUE``)`` `` ``# collect two samples`` ``sim``$``sample_cells``(``"S1"``, `[`c`](https://rdrr.io/r/base/c.html)`(``145``, ``230``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``215``, ``300``)``)`` ``sim``$``sample_cells``(``"S2"``, `[`c`](https://rdrr.io/r/base/c.html)`(``350``, ``300``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``420``, ``370``)``)`` `` ``# add a further mutant and derive it from B`` ``sim``$``add_mutant``(``name ``=`` ``"D"``, `[`c`](https://rdrr.io/r/base/c.html)`(``duplication ``=`` ``0.8``, death ``=`` ``0.01``)``)`` ``sim``$``mutate_progeny``(``sim``$``choose_border_cell_in``(``"B"``)``, ``"D"``)`` `` ``# let the tumour evolve until the mutant C and D cumulatively`` ``# consist of 10000 cells`` ``sim``$``run_until``(``sim``$``var``(``"C"``)`` ``+`` ``sim``$``var``(``"D"``)`` ``==`` ``1e5``, quiet ``=`` ``TRUE``)`` `` ``# collect the last sample`` ``sim``$``sample_cells``(``"S3"``, `[`c`](https://rdrr.io/r/base/c.html)`(``100``, ``400``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``170``, ``470``)``)`` `` ``# build the forest`` ``sample_forest`` ``<-`` ``sim``$``get_sample_forests``(``)`

At the end of the evolution, by calling

[`build_snapshot_video`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)`(``sim``, res ``=`` ``300``, width``=``1440``, height``=``1024``, framerate ``=`` ``20``)`

we produce the following video.

## Customizing Pauses

Users can specify pauses on events declaring the length of the pause on
specific sample or new mutant. They can also add pauses on any sample or
new mutant.

``` r
pauses_on_event <- list(
    "mutant emerged" = list(

        # pause 4 seconds when "C" arises
        "C" = as.difftime(4, units="secs"),

        # pause 4 frames when "B" arises
        "B" = 4),

    # pause 3 seconds at each sample
    "sampling" = as.difftime(3, units="secs")
)

build_snapshot_video(sim, res = 300, width=1440, height=1024, framerate = 20
                     pauses_on_event = pauses_on_event)
```

Pauses can also applied to specific frames.

``` r
pauses_on_frame <- list(
    # pause 3 seconds on frame 101
    list(frame = 101, length = as.difftime(3, units="secs")),

    # pause 3 seconds on frame 102
    list(frame = 102, length = as.difftime(3, units="secs"))
)

build_snapshot_video(sim, res = 300, width=1440, height=1024, framerate = 20
                     pauses_on_frame = pauses_on_frame)
```

## Customizing Frames

[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)
also supports frame customisations. The parameter `frame_generator`
allows users to specify their frame plotting function. The actual
parameter should be a function that takes in input a simulation snapshot
and returns a plot corresponding to a frame of the video.

For instance, we can produce a video focusing on the cells in a
rectangle. First of all, we need to define a focus function for the
cells in the rectangle (see,
[`plot_tissue()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)),
and, then, we changes define the `frame_generator()`.

`# a focus function that highlights the cells in the rectangle`` ``# [250,350]x[300,350]`` ``in_rectangle`` ``<-`` ``function``(``cells``)`` ``{`` `` ``(``cells``$``position_x`` ``>=`` ``250`` ``&`` ``cells``$``position_x`` ``<=`` ``350`` `` ``&`` ``cells``$``position_y`` ``>=`` ``300`` ``&`` ``cells``$``position_y`` ``<=`` ``350``)`` ``}`` `` ``# define a unique color map for all the frames`` ``color_map`` ``<-`` `[`get_species_colors`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_species_colors.md)`(``sim``)`` `` ``# this function will be called to plot each video frames`` ``plot_in_rectangle`` ``<-`` ``function``(``snapshot``)`` ``{`` `` `[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``snapshot``,`` `` color_map ``=`` ``color_map``,`` `` plot_sample_region ``=`` ``FALSE``,`` `` list_all_labels ``=`` ``TRUE``,`` `` focus_function ``=`` ``in_rectangle``)`` ``}`

The snippet

[`build_snapshot_video`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)`(``sim``, res ``=`` ``300``, width``=``1440``, height``=``1024``,`` `` framerate ``=`` ``20``, frame_generator ``=`` ``plot_in_rectangle``)`

produces the following result.

## Predefined Frame Generators

`ProCESS` provides some default frame generator using the functions
\[get_tissue_frame_gen()\], \[get_forest_frame_gen()\], and
\[get_tissue_forest_frame_gen()\]. They all produces a named list of two
functions: `frame_generator` and `cleanup_function`. The former is a
frame generator; the latter a function that removes the temporary files,
and that must be called after the frame generation process has been
completed. The function \[get_tissue_frame_gen()\] produces a frame
generator that represents the tissue configuration of a snapshot in each
frame. This frame generator is the default one for
\[build_snapshot_video()\]. The function \[get_forest_frame_gen()\]
returns a frame generator depicting each snapshot as the sample forest
at the snapshot time. It requires a temporary file and the cleanup
function returned by \[get_forest_frame_gen()\] removes it. Finally, the
function \[get_forest_frame_gen()\] returns a frame generator that
represents each snapshot in the simulation as both its tissue
configuration and the sample forest at the snapshot time. This function
is only available when \[patchwork::patchwork-package\] is installed in
the system.

## Parallel Frame Generation

When the packages [furrr](https://CRAN.R-project.org/package=furrr) and
[progressr](https://CRAN.R-project.org/package=progressr) are installed
in the system, the function
[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)
uses them to parallelise frame generations. However, ProCESS data
structures are not serializable in R and cannot be transparently passed
to the [furrr](https://CRAN.R-project.org/package=furrr) processes.
Because of that, all ProCESS objects involved in the evaluation of
`frame_generator()`, with the exception of the function parameter, must
be saved to the disk outside `frame_generator()` and, then, reloaded
inside it. In this way, the objects will be local to each
[furrr](https://CRAN.R-project.org/package=furrr) process.

For instance, if we want to generate a video that focus on the cells
represented in the sample forest, we first need to save the sample
forest.

`sample_forest``$``save``(``"forest.sff"``)`` ``#> [█---------------------------------------] 0% [00m:00s] Saving forest [████████████████████████████████████████] 100% [00m:00s] Forest saved`

Then, we must define a `frame_generator()` that loads its own copy of
the sample forest and use it.

`# this function will be called to plot each video frames`` ``plot_in_forest`` ``<-`` ``function``(``snapshot``)`` ``{`` `` ``# load the local copy of the forest`` `` ``local_forest`` ``<-`` `[`load_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_forest.md)`(``"forest.sff"``, quiet ``=`` ``TRUE``)`` `` `` ``# define a unique color map for all the frames`` `` ``color_map`` ``<-`` `[`get_species_colors`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_species_colors.md)`(``local_forest``)`` `` `` ``# a focus function that highlights the cells represented by`` `` ``# the sample forest`` `` ``in_sample_forest`` ``<-`` ``function``(``cells``)`` ``{`` `` ``# sample_forest has been defined in four_mutants.md`` `` ``local_forest``$``represents_cell``(``cells``$``cell_id``)`` `` ``}`` `` `` ``# plot the tissue using the shared color map `` `` `[`plot_tissue`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)`(``snapshot``,`` `` color_map ``=`` ``color_map``,`` `` plot_sample_region ``=`` ``FALSE``,`` `` list_all_labels ``=`` ``TRUE``,`` `` focus_function ``=`` ``in_sample_forest``)`` ``}`

Please, notice that `color_map` is a serialisable object and is
automatically passed to all frame generator processes.

Finally, the function call

[`build_snapshot_video`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)`(``sim``, res ``=`` ``300``, width``=``1440``, height``=``1024``,`` `` framerate ``=`` ``20``, frame_generator ``=`` ``plot_in_forest``,`` `` pauses_on_event ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``)``)`

produces the following video.

The package [patchwork](https://CRAN.R-project.org/package=patchwork)
can be used to generate more complex videos.

The `frame_generator` built by \[get_tissue_forest_frame_gen()\]
exploits this package.

`fg`` ``<-`` `[`get_tissue_forest_frame_gen`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_tissue_forest_frame_gen.md)`(``sim``)`` `` `[`build_snapshot_video`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)`(``sim``, width ``=`` ``1820``, height ``=`` ``1024``, framerate ``=`` ``20``,`` `` res ``=`` ``300``, frame_generator ``=`` ``fg``$``frame_generator``)`` `` ``fg``$``cleanup_function``(``)`
