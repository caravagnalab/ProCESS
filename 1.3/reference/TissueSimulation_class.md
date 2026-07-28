# Simulating the cell evolution in a tissue

This class simulates the cell evolution on a tissue.

## Details

The objects of this class can simulate the evolution of many cells
belonging to different *species* on a tissue. Each cell can duplicate or
die according to the rates that delineate the cell species.

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)
supports epigenetic evolutions, and it lets users define species pairs
that belong to the same mutant (even though, its genomic
characterisation is unknown) and differ because of their epigenetic
states.

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)
models epigenetic mutations and allows a cell in one of mutant species
to generate a new cell belonging to the other species of the same mutant
at a specified rate.

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)
also allows users to schedule mutations from one mutant to a different
mutant.

This class provides the following methods and properties:

- [`add_epigenetic_state()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_state.md)
  adds an epigenetic state and its species.

- [`add_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_states.md)
  adds epigenetic states and their species.

- [`add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutant.md)
  adds a mutant and its species.

- [`add_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutants.md)
  adds mutants and their species.

- [`border_growth_model`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-border_growth_model.md)
  is a Boolean flag to enable/disable the border-driven growth model.

- [`choose_border_cell_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-choose_border_cell_in.md)
  randomly chooses a cell in tumour border.

- [`choose_cell_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-choose_cell_in.md)
  randomly chooses a tumour cell.

- [`death_activation_level`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-death_activation_level.md)
  stores the minimum number of cells in a species to enable death.

- [`get_added_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_added_cells.md)
  returns the cells manually added to the simulation.

- [`get_cell()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_cell.md)
  returns one of the simulation cells by its position.

- [`get_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_cells.md)
  returns some of the simulation cells by their positions.

- [`get_clock()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_clock.md)
  returns the current simulation time.

- [`get_count_history()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_count_history.md)
  returns the number of cells by species along the simulation.

- [`get_counts()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_counts.md)
  returns the number of cells by species.

- [`get_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_epigenetic_states.md)
  returns the names of the simulated epigenetic states.

- [`get_firing_history()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_firing_history.md)
  returns the number of events fired along the simulation.

- [`get_firings()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_firings.md)
  returns the number of fired events.

- [`get_lineage_graph()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_lineage_graph.md)
  returns the lineage graph.

- [`get_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_mutants.md)
  returns the names of the simulated mutants.

- [`get_name()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_name.md)
  returns the name of the simulation.

- [`get_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_rates.md)
  returns the tissue simulation rates.

- [`get_rates_update_history()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_rates_update_history.md)
  returns the history of the rates along the simulation.

- [`get_sample_forest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_sample_forest.md)
  returns the sample forest of the tissue simulation.

- [`get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_samples_info.md)
  returns information about the collected samples.

- [`get_snapshot_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_snapshot_info.md)
  returns information about the collected snapshots.

- [`get_tissue_size()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_tissue_size.md)
  returns the size of the simulated tissue.

- [`history_delta`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-history_delta.md)
  stores the time difference between two consecutive samples of the
  tissue simulation status.

- [`mutate_progeny()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-mutate_progeny.md)
  duplicates a cell and mutates one of its children.

- [`place_cell()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-place_cell.md)
  places a cell in the simulated tissue.

- [`run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)
  simulates the evolution of the tissue until a condition holds.

- [`run_up_to_event()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_event.md)
  simulates the evolution of the tissue until the number of events of
  either a mutant or a species is below a specified threshold.

- [`run_up_to_size()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_size.md)
  simulates the evolution of the tissue until the number of cells of
  either a mutant or a species is below a specified threshold .

- [`run_up_to_time()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_up_to_time.md)
  simulates the evolution of the tissue until the specified simulation
  time.

- [`sample_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-sample_cells.md)
  samples tissue cells.

- [`schedule_mutation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-schedule_mutation.md)
  schedules a mutation.

- [`search_sample()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-search_sample.md)
  searches for a rectangular tissue sample satisfying some conditions.

- [`search_samples()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-search_samples.md)
  searches for rectangular tissue samples satisfying some conditions.

- [`set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md)
  sets the rate of an event.

- [`set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md)
  sets the rates of multiple events.

- [`snapshot_trigger`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_trigger.md)
  exposes the tissue simulation snapshot triggers.

- [`var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md)
  builds a simulation variable.

## See also

[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)
