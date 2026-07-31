/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2026 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>

#include <Rcpp.h>

#include "sample_forest.hpp"
#include "simulation.hpp"
#include "tissue_rectangle.hpp"

#include "node_tour.hpp"

using namespace Rcpp;

namespace RE = CLONES::Mutants::Evolutions;
namespace RC = CLONES::Mutants;

RCPP_MODULE(Mutants)
{

//' @name TissueRectangle
//' @title A rectangle in the tissue
//' @details The objects of this class represent a rectangle of tissue.
//'   They provides the following properties and methods:
//'   - <code>[TissueRectangle$lower_corner]</code> stores the lower
//'     left corner of the rectangle.
//'   - <code>[TissueRectangle$upper_corner]</code> stores the upper
//'     right corner of the rectangle.
    class_<TissueRectangle>("TissueRectangle")
//' @name TissueRectangle$new
//' @title Build a new rectangle of tissue
//' @examples
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), c(550, 475))
//'
//' rect
//'
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), 50, 25)
//'
//' rect
        .constructor<std::vector<uint16_t>, std::vector<uint16_t>>(
            "Create a new rectangle")
        .constructor<std::vector<uint16_t>, RE::AxisSize, RE::AxisSize>(
            "Create a new rectangle")

//' @name TissueRectangle$lower_corner
//' @title The tissue rectangle lower corner
//' @description This is the vector of the minima among all the rectangle
//'   dimensions.
//' @examples
//' # create a 51x51-rectangle from (500, 500) to (550, 550)
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the rectangle lower corner, i.e., (500, 500)
//' rect$lower_corner
//' @seealso <code>[TissueRectangle]</code>
        .property("lower_corner", &TissueRectangle::get_lower_corner,
                  "The rectangle lower corner")

//' @name TissueRectangle$upper_corner
//' @title The tissue rectangle upper corner
//' @description This is the vector of the maxima among all the rectangle
//'   dimensions.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the rectangle upper corner, i.e., (550, 550)
//' rect$upper_corner
//' @seealso <code>[TissueRectangle]</code>
        .property("upper_corner", &TissueRectangle::get_upper_corner,
                  "The rectangle upper corner")
        .method("show", &TissueRectangle::show);

//' @name TissueSimulation_class
//' @title Simulating the cell evolution in a tissue
//' @description This class simulates the cell evolution on a tissue.
//' @details The objects of this class can simulate the evolution
//'   of many cells belonging to different *species* on a tissue. Each
//'   cell can duplicate or die according to the rates that delineate
//'   the cell species.
//'
//'   [TissueSimulation()] supports epigenetic evolutions, and it lets
//'   users define species pairs that belong to the same mutant (even
//'   though, its genomic characterisation is unknown) and differ
//'   because of their epigenetic states.
//'
//'   [TissueSimulation()] models epigenetic mutations and allows a cell
//'   in one of mutant species to generate a new cell belonging to the
//'   other species of the same mutant at a specified rate.
//'
//'   [TissueSimulation()] also allows users to schedule mutations from
//'   one mutant to a different mutant.
//'
//'   This class provides the following methods and properties:
//'   - <code>[add_epigenetic_state()](TissueSimulation-cash-add_epigenetic_state.md)</code>
//'     adds an epigenetic state and its species.
//'   - <code>[add_epigenetic_states()](TissueSimulation-cash-add_epigenetic_states.md)</code>
//'     adds epigenetic states and their species.
//'   - <code>[add_mutant()](TissueSimulation-cash-add_mutant.md)</code>
//'     adds a mutant and its species.
//'   - <code>[add_mutants()](TissueSimulation-cash-add_mutants.md)</code>
//'     adds mutants and their species.
//'   - <code>[border_growth_model](TissueSimulation-cash-border_growth_model.md)</code>
//'     is a Boolean flag to enable/disable the border-driven growth model.
//'   - <code>[choose_border_cell_in()](TissueSimulation-cash-choose_border_cell_in.md)</code>
//'     randomly chooses a cell in tumour border.
//'   - <code>[choose_cell_in()](TissueSimulation-cash-choose_cell_in.md)</code>
//'     randomly chooses a tumour cell.
//'   - <code>[death_activation_level](TissueSimulation-cash-death_activation_level.md)</code>
//'     stores the minimum number of cells in a species to enable death.
//'   - <code>[get_added_cells()](TissueSimulation-cash-get_added_cells.md)</code>
//'     returns the cells manually added to the simulation.
//'   - <code>[get_cell()](TissueSimulation-cash-get_cell.md)</code>
//'     returns one of the simulation cells by its position.
//'   - <code>[get_cells()](TissueSimulation-cash-get_cells.md)</code>
//'     returns some of the simulation cells by their positions.
//'   - <code>[get_clock()](TissueSimulation-cash-get_clock.md)</code>
//'     returns the current simulation time.
//'   - <code>[get_count_history()](TissueSimulation-cash-get_count_history.md)</code>
//'     returns the number of cells by species along the simulation.
//'   - <code>[get_counts()](TissueSimulation-cash-get_counts.md)</code>
//'     returns the number of cells by species.
//'   - <code>[get_epigenetic_states()](TissueSimulation-cash-get_epigenetic_states.md)</code>
//'     returns the names of the simulated epigenetic states.
//'   - <code>[get_firing_history()](TissueSimulation-cash-get_firing_history.md)</code>
//'     returns the number of events fired along the simulation.
//'   - <code>[get_firings()](TissueSimulation-cash-get_firings.md)</code>
//'     returns the number of fired events.
//'   - <code>[get_lineage_graph()](TissueSimulation-cash-get_lineage_graph.md)</code>
//'     returns the lineage graph.
//'   - <code>[get_mutants()](TissueSimulation-cash-get_mutants.md)</code>
//'     returns the names of the simulated mutants.
//'   - <code>[get_name()](TissueSimulation-cash-get_name.md)</code>
//'     returns the name of the simulation.
//'   - <code>[get_rates()](TissueSimulation-cash-get_rates.md)</code>
//'     returns the tissue simulation rates.
//'   - <code>[get_rates_update_history()](TissueSimulation-cash-get_rates_update_history.md)</code>
//'     returns the history of the rates along the simulation.
//'   - <code>[get_sample_forest()](TissueSimulation-cash-get_sample_forest.md)</code>
//'     returns the sample forest of the tissue simulation.
//'   - <code>[get_samples_info()](TissueSimulation-cash-get_samples_info.md)</code>
//'     returns information about the collected samples.
//'   - <code>[get_snapshot_info()](TissueSimulation-cash-get_snapshot_info.md)</code>
//'     returns information about the collected snapshots.
//'   - <code>[get_tissue_size()](TissueSimulation-cash-get_tissue_size.md)</code>
//'     returns the size of the simulated tissue.
//'   - <code>[history_delta](TissueSimulation-cash-history_delta.md)</code>
//'     stores the time difference between two consecutive samples
//'     of the tissue simulation status.
//'   - <code>[mutate_progeny()](TissueSimulation-cash-mutate_progeny.md)</code>
//'     duplicates a cell and mutates one of its children.
//'   - <code>[place_cell()](TissueSimulation-cash-place_cell.md)</code>
//'     places a cell in the simulated tissue.
//'   - <code>[run_until()](TissueSimulation-cash-run_until.md)</code>
//'     simulates the evolution of the tissue until a condition holds.
//'   - <code>[run_up_to_event()](TissueSimulation-cash-run_up_to_event.md)</code>
//'     simulates the evolution of the tissue until the number of events
//'     of either a mutant or a species is below a specified threshold.
//'   - <code>[run_up_to_size()](TissueSimulation-cash-run_up_to_size.md)</code>
//'     simulates the evolution of the tissue until the number of cells
//'     of either a mutant or a species is below a specified threshold .
//'   - <code>[run_up_to_time()](TissueSimulation-cash-run_up_to_time.md)</code>
//'     simulates the evolution of the tissue until the specified simulation
//'     time.
//'   - <code>[sample_cells()](TissueSimulation-cash-sample_cells.md)</code>
//'     samples tissue cells.
//'   - <code>[schedule_mutation()](TissueSimulation-cash-schedule_mutation.md)</code>
//'     schedules a mutation.
//'   - <code>[search_sample()](TissueSimulation-cash-search_sample.md)</code>
//'     searches for a rectangular tissue sample satisfying some conditions.
//'   - <code>[search_samples()](TissueSimulation-cash-search_samples.md)</code>
//'     searches for rectangular tissue samples satisfying some conditions.
//'   - <code>[set_rate()](TissueSimulation-cash-set_rate.md)</code>
//'     sets the rate of an event.
//'   - <code>[set_rates()](TissueSimulation-cash-set_rates.md)</code>
//'     sets the rates of multiple events.
//'   - <code>[snapshot_trigger](TissueSimulation-cash-snapshot_trigger.md)</code>
//'     exposes the tissue simulation snapshot triggers.
//'   - <code>[var()](TissueSimulation-cash-var.md)</code>
//'     builds a simulation variable.
//' @keywords internal
//' @seealso [TissueSimulation()]
    class_<TissueSimulation>("TissueSimulation")

//' @name TissueSimulation$place_cell
//' @title Placing one cell in the tissue
//' @description This method places a cell in the tissue.
//' @param species The name of the new cell species.
//' @param x The position on the x axis of the cell.
//' @param y The position on the y axis of the cell.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
//'
//' # place a cell of species "A" in position (500,500)
//' sim$place_cell("A", 500, 500)
//' @examples
//' # create a simulation
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.1)))
//'
//' # place a cell of species "A[E1]" in position (500,500)
//' sim$place_cell("A[E1]", 500, 500)
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("place_cell", (void (TissueSimulation::*)(const SEXP&, const SEXP&, const SEXP&))
                &TissueSimulation::place_cell, "Placing a cell in the tissue")

//' @name TissueSimulation$add_mutant
//' @title Adding a mutant and its species
//' @description This method adds a mutant and its species to the
//'   simulation.
//' @details This method adds a mutant to the simulation. The method also
//'   creates the species associated to the new mutant according to the known
//'   epigenetic states. The default rate of the new species is set to zero.
//'   Optionally, user can provide a list specifying the rates of the
//'   associated species.
//' @param mutant The mutant name.
//' @param rate_list The list of the mutant's rates (optional).
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # see the simulation setup
//' sim
//'
//' # add the mutant "A" to the simulation.
//' sim$add_mutant(name = "A")
//'
//' # see the simulation setup
//' sim
//'
//' # add the mutant "B" to the simulation and set its duplication and death rates
//' sim$add_mutant(name = "B", rate_list = c(duplication = 0.3, death = 0.1))
//'
//' # see the simulation setup
//' sim
//'
//' # add epigenetic states (rates are reset)
//' sim$add_epigenetic_states(c("E1", "E2", "E3"))
//'
//' # add the mutant "C" to the simulation, set the duplication and death rates
//' # of all its species, and differentiate "C[E1]" by setting its death rate
//' # and the rates of the switch toward "C[E2]" and "C[E3]".
//' sim$add_mutant("C", list(duplication = 0.3, death = 0.1,
//'                          E1=list(death = 0.2, E2=0.01, E3=0.1)))
//'
//' # see the simulation setup
//' sim
//' @seealso [TissueSimulation$add_mutants()], [TissueSimulation$add_epigenetic_state()],
//'   [TissueSimulation$add_epigenetic_states()], [TissueSimulation$get_mutants()],
//'   [TissueSimulation$set_rate()], [TissueSimulation$set_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "add_mutant",
            (void (TissueSimulation::*)(const SEXP &))(&TissueSimulation::add_mutant),
            "Add a new mutant")
        .method(
            "add_mutant",
            (void (TissueSimulation::*)(const SEXP &, const SEXP&))(&TissueSimulation::add_mutant),
            "Add a new mutant")

//' @name TissueSimulation$add_mutants
//' @title Adding mutants and their species
//' @description This method adds mutants and their species to the
//'   simulation.
//' @details This method adds mutants to the simulation. The method also
//'   creates the species associated to the new mutants according to the known
//'   epigenetic states. The default rate of the new species is set to zero.
//' @param mutants A list of mutant names
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' sim$get_mutants()
//'
//' # add the mutants "A", "B", and "C" to the simulation.
//' sim$add_mutants(c("A", "B", "C"))
//'
//' sim$get_mutants()
//' @seealso [TissueSimulation$add_mutant()], [TissueSimulation$add_epigenetic_state()],
//'   [TissueSimulation$add_epigenetic_states()], [TissueSimulation$get_mutants()],
//'   [TissueSimulation$set_rate()], [TissueSimulation$set_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "add_mutants",
            (void (TissueSimulation::*)(const std::list<std::string> &))(&TissueSimulation::add_mutants),
            "Add new mutants")

//' @name TissueSimulation$add_epigenetic_state
//' @title Adding an epigenetic state and its species
//' @description This method adds an epigenetic state and its species to the
//'   simulation.
//' @details This method introduces a new epigenetic state into the
//'   simulation. Additionally, the method adds to each known mutant a new
//'   species. The default rate of the new species is set to zero.
//' @param epigenetic_state The name of the epigenetic state to add.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' sim$get_epigenetic_states()
//'
//' # add the epigenetic state "E1" to the simulation.
//' sim$add_epigenetic_state("E1")
//'
//' sim$get_epigenetic_states()
//' @seealso [TissueSimulation$add_epigenetic_states()],
//'   [TissueSimulation$get_epigenetic_states()],
//'   [TissueSimulation$add_mutant()], [TissueSimulation$add_mutants()],
//'   [TissueSimulation$set_rate()], [TissueSimulation$set_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "add_epigenetic_state",
            (void (TissueSimulation::*)(const std::string &))(&TissueSimulation::add_epistate),
            "Add a new epigenetic space")

//' @name TissueSimulation$add_epigenetic_states
//' @title Adding epigenetic states and their species
//' @description This method adds epigenetic states and their species to the
//'   simulation.
//' @details This method introduces novel epigenetic states into the
//'   simulation. Additionally, the method adds to each known mutant as many
//'   species as the number of newly introduced epigenetic states. The default
//'   rate of these new species is set to zero.
//' @param epigenetic_states A list of epigenetic state names
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' sim$get_epigenetic_states()
//'
//' # add the epigenetic state "E1", "E2", and "E3" to the simulation.
//' sim$add_epigenetic_states(c("E1", "E2", "E3"))
//'
//' sim$get_epigenetic_states()
//' @seealso [TissueSimulation$add_epigenetic_state()],
//'   [TissueSimulation$get_epigenetic_states()],
//'   [TissueSimulation$add_mutant()], [TissueSimulation$add_mutants()],
//'   [TissueSimulation$set_rate()], [TissueSimulation$set_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "add_epigenetic_states",
            (void (TissueSimulation::*)(const std::list<std::string> &))(&TissueSimulation::add_epistates),
            "Add new epigenetic spaces")

//' @name TissueSimulation$choose_cell_in
//' @title Picking one border cell
//' @description This method chooses one cell among those
//'   belonging to one of the specified mutants and species.
//' @details It randomly chooses one of the cells belonging to
//'   either mutants or species. Optionally, the lower and upper
//'   corners of a tissue rectangular selection can be provided to
//'   obtain one cell in the rectangle.
//' @param names The names of the mutants or species among which the
//'   cell must be choosen. Can either be a single name or a list of
//'   names.
//' @param lower_corner The lower corner of the rectangular selection
//'   (optional).
//' @param upper_corner The upper corner of the rectangular selection
//'   (optional).
//' @return A list reporting `cell_id`, `mutant`, `position_x`,
//'   `position_y`, and, when the simulation has epigenetic states,
//'   `epistate` of the chosen cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.1)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.15, death = 0.1, E2 = 0.1),
//'                     E2 = list(duplication = 0.4, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 20)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # set the death activation level
//' sim$death_activation_level <- 100
//'
//' # run the simulation until "B[E2]" accounts for less than 1000 cells
//' sim$run_up_to_size("B[E2]", 1000)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_cell_in("B")
//'
//' # Randomly choose one cell in "B" in a rectangle
//' sim$choose_cell_in("B", c(500, 500), c(520, 520))
//'
//' # Randomly choose one cell in "B[E1]"
//' sim$choose_cell_in("B[E1]")
//'
//' # Randomly choose one cell in "B[E1]"
//' sim$choose_cell_in(c("B[E1]", "A"))
//' @seealso [TissueSimulation$choose_border_cell_in()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("choose_cell_in",
                (List (TissueSimulation::*)(const SEXP &))(
                    &TissueSimulation::choose_cell_in),
                "Randomly choose one cell in a mutant")
        .method("choose_cell_in",
                (List (TissueSimulation::*)(const SEXP &,
                                            const std::vector<RE::AxisPosition> &,
                                            const std::vector<RE::AxisPosition> &))(
                    &TissueSimulation::choose_cell_in),
                "Randomly choose one cell having a specified mutant in a rectangular "
                "selection")

//' @name TissueSimulation$choose_border_cell_in
//' @title Picking one border cell
//' @description This method chooses one border cell among those
//'   belonging to one of the specified mutants and species.
//' @details It randomly chooses one of the cells belonging to
//'   either a mutant or a species that has a wild-type cell in its
//'   neighborhood. Optionally, the lower and upper corners of a
//'   tissue rectangular selection can be provided to obtain one
//'   cell in the rectangle.
//' @param names The names of the mutants or species among which the
//'   cell must be choosen. Can either be a single name or a list of
//'   names.
//' @param lower_corner The lower corner of the rectangular selection
//'   (optional).
//' @param upper_corner The upper corner of the rectangular selection
//'   (optional).
//' @return A list reporting `cell_id`, `mutant`, `position_x`,
//'   `position_y`, and, when the simulation has epigenetic states,
//'   `epistate` of the chosen cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.1)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.15, death = 0.1, E2 = 0.1),
//'                     E2 = list(duplication = 0.4, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 20)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # set the death activation level
//' sim$death_activation_level <- 100
//'
//' # run the simulation until "B[E2]" accounts for less than 1000 cells
//' sim$run_up_to_size("B[E2]", 1000)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_border_cell_in("B")
//'
//' # Randomly choose one cell in "B" in a rectangle
//' sim$choose_border_cell_in("B", c(500, 500), c(520, 520))
//'
//' # Randomly choose one cell in "B[E1]"
//' sim$choose_border_cell_in("B[E1]")
//'
//' # Randomly choose one cell in "B[E1]" and any species in "A"
//' sim$choose_border_cell_in(c("B[E1]", "A"))
//' @seealso [TissueSimulation$choose_cell_in()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("choose_border_cell_in",
                (List (TissueSimulation::*)(const SEXP &))(
                    &TissueSimulation::choose_border_cell_in),
                "Randomly choose one cell in a mutant")
        .method("choose_border_cell_in",
                (List (TissueSimulation::*)(const SEXP &,
                                            const std::vector<RE::AxisPosition> &,
                                            const std::vector<RE::AxisPosition> &))(
                    &TissueSimulation::choose_border_cell_in),
                "Randomly choose one cell having a specified mutant in a rectangular "
                "selection")

//' @name TissueSimulation$snapshot_triggers
//' @title The snapshot triggers
//' @description This property is the tissue simulation snapshot trigger list
//' @details This property is a named list containing three values at most:
//'   `time interval`, `clock interval`, and `number of cells`. They
//'   represents the maximum computation time, the maximum simulation time,
//'   and the maximum difference in the number of tumour cells between two
//'   snapshots, respectively.
//'
//'   Notice that this property differs from
//'   <code>[TissueSimulation$history_delta](TissueSimulation-cash-history_delta.md)</code>.
//' @examples
//' # set the random seed
//' set.seed(0)
//'
//' # build a simulation
//' sim <- TissueSimulation()
//'
//' # add a mutant
//' sim$add_mutant("A", list(duplication = 3, death = 1))
//'
//' # place a cell
//' sim$place_cell("A", 500, 500)
//'
//' # get snapshot information
//' sim$get_snapshot_info()
//'
//' # get the simulation's snapshot triggers
//' sim$snapshot_triggers
//'
//' # take a new snapshot every 1000 new tumour cells
//' sim$snapshot_triggers <- list("number of cells" = 1000)
//'
//' # get new simulation's snapshot triggers
//' sim$snapshot_triggers
//'
//' # let the simulation evolve until consists of 5000 cells
//' sim$run_up_to_size("A", 4000)
//'
//' # get snapshot information
//' sim$get_snapshot_info()
//'
//' # set new snapshot triggers: get a snapshot every 30 seconds, every 10
//' # simulated time units, and every 1000 new cells (because set previously)
//' sim$snapshot_triggers <- list("time interval" = as.difftime(30, units = "secs"),
//'                               "clock interval" = 10)
//'
//' # get current simulation's snapshot triggers
//' sim$snapshot_triggers
//'
//' # take a snapshot every 5 simulated time units and reset the other
//' # trigger types
//' sim$snapshot_triggers <- list("time interval" = NULL,
//'                               "clock interval" = 5,
//'                               "number of cells" = NULL)
//'
//' # get new simulation's snapshot triggers
//' sim$snapshot_triggers
//'
//' # let the simulation evolve for other 20 time units
//' sim$run_up_to_time(sim$get_clock() + 20)
//'
//' # get snapshot information
//' sim$get_snapshot_info()
//' @seealso <code>[TissueSimulation](TissueSimulation_code.md)</code>,
//'   [TissueSimulation$get_snapshot_info()],
//'   [recover_simulation()],
//'   <code>[TissueSimulation$history_delta](TissueSimulation-cash-history_delta.md)</code>
        .property("snapshot_triggers",
                  &TissueSimulation::get_snapshot_triggers,
                  &TissueSimulation::set_snapshot_triggers,
                 "The tissue simulation snapshot triggers")

//' @name TissueSimulation$get_snapshot_info
//' @title Getting snapshot information
//' @description This method returns the data frame of the snapshots
//' @return A data frame consisting of four columns: `time`, `clock`,
//'   `cells`, and `file`. Each row represents a snapshot. The column
//'   `time` stores the snapshot time. The columns `clock` and `cells`
//'   maintain the simulated time and the number of tumour cells at
//'   the snapshot time, respectively. Finally, the column `file`
//'   contains the snapshot file path.
//' @examples
//' # set the random seed
//' set.seed(0)
//'
//' # build a simulation
//' sim <- TissueSimulation()
//'
//' # add a mutant
//' sim$add_mutant("A", list(duplication = 3, death = 1))
//'
//' # place a cell
//' sim$place_cell("A", 500, 500)
//'
//' # get snapshot information
//' sim$get_snapshot_info()
//'
//' # take a new snapshot every 10 simulated time units
//' sim$snapshot_triggers <- list("clock interval" = 10)
//'
//' # let the simulation evolve up to time 30
//' sim$run_up_to_time(30)
//'
//' # get snapshot information
//' sim$get_snapshot_info()
//' @seealso <code>[TissueSimulation](TissueSimulation_code.md)</code>,
//'   <code>[TissueSimulation$snapshot_triggers](TissueSimulation-cash-snapshot_triggers.md)</code>,
//'   [recover_simulation()]
        .method("get_snapshot_info",
                &TissueSimulation::get_snapshot_dataframe,
                "Get a data frame of the snapshots")

//' @name TissueSimulation$schedule_mutation
//' @title Scheduling a mutation
//' @description This method schedules a mutant mutation
//' @details The mutation can occur from any of the species of
//'   the source mutant to the species of the destination mutant
//'   with a consistent epigenetic state.
//'   For the sake of example, if the simulation has no epigenetic
//'   states and a mutation from "A" to "B" is scheduled, then
//'   the first duplication of an "A"'s cell occurring after the
//'   specified time generates two cells: one of them belong to "A"
//'   and the other to "B".
//'   Analogously, if the simulation has epigenetic states and the
//'   first duplication of an "A"'s cell after the specified time
//'   occurs to a cell in the species "A[Ei]", then the offspring
//'   of the cell consists of one cell in "A[Ei]" and one in "B[Ei]".
//' @param src The name of the mutant from which the mutation occurs.
//' @param dest The name of the mutant to which the mutation leads.
//' @param time The simulated time at which the mutation will occurs.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", c(duplication = 0.4, death = 0.01))
//'
//' # schedule an evolution from mutant "A" to mutant "B" at time 50
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//'
//' # place a cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to the first cell in B
//' sim$run_up_to_size("B", 1)
//'
//' sim
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("schedule_mutation", &TissueSimulation::schedule_mutation,
                "Add a timed mutation between two different species")

//' @name TissueSimulation$get_mutants
//' @title Getting the names of the simulated mutants
//' @description This method returns the names of the simulated mutants.
//' @return A data frame containing the column `mutant`. Each row of the data frame reports the
//'   name of one of the simulated mutants.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' sim$get_mutants()
//'
//' # add the mutants "A", "B", and "C" to the simulation.
//' sim$add_mutants(c("A", "B", "C"))
//'
//' sim$get_mutants()
//' @seealso [TissueSimulation$add_mutant()], [TissueSimulation$add_mutants()],
//'   [TissueSimulation$get_epigenetic_states()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_mutants", &TissueSimulation::get_mutant_names,
                "Get the mutants added to the simulation")

//' @name TissueSimulation$get_epigenetic_states
//' @title Getting the epigenetic states in the simulation
//' @description This method returns the epigenetic states in the simulations.
//' @return A data frame having a single column `epistate`. The column
//'   contains the names of the epigenetic states added to the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # get the data frame of the mutants
//' sim$get_epigenetic_states()
//'
//' # add epigenetic states
//' sim$add_epigenetic_states(c("E1","E2"))
//'
//' # get the data frame of the mutants
//' sim$get_epigenetic_states()
//' @seealso [TissueSimulation$add_epigenetic_state()],
//'   [TissueSimulation$add_epigenetic_states()]
//'   [TissueSimulation$get_mutants()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_epigenetic_states", &TissueSimulation::get_epigenetic_state_names,
                "Get the epigenetic states added to the simulation")

//' @name TissueSimulation$get_sample_forest
//' @title Getting the sample forest
//' @description This method returns the sample forest.
//' @return The sample forest having as leaves the sampled cells
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size("A", 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", c(450,475), c(500,550))
//'
//' # build the sample forest
//' forest <- sim$get_sample_forest()
//'
//' forest
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_sample_forest", &TissueSimulation::get_sample_forest,
                "Get the sample forest having as leaves the sampled cells")

//' @name TissueSimulation$death_activation_level
//' @title Death activation level
//' @description The number of cells that activates cell death in a species.
//' @details This value is the minimum number of cells that
//'   enables cell death in a species. The cell of a species $S$ can die
//'   if and only if that $S$ has reached the death activation level at
//'   least once during the simulation.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # get the simulation death activation level
//' sim$death_activation_level
//'
//' # set the death activation level to 50
//' sim$death_activation_level <- 50
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .property("death_activation_level", &TissueSimulation::get_death_activation_level,
                  &TissueSimulation::set_death_activation_level,
                  "The number of cells in a species that activates cell death")

//' @name TissueSimulation$border_growth_model
//' @title Internal cells duplication
//' @description This property switches between homogeneous and border driven
//'   growth models.
//' @details This Boolean flag switches between homogeneous and border driven
//'   growth models. When it is set to `TRUE`, the border-growth model
//'   is used. Otherwise, the homogeneous-growth model is applied.
//'   It is set to `TRUE` by default.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # is the simulation using the border driven growth model
//' # (default: `TRUE`)
//' sim$border_growth_model
//'
//' # switch to the homogeneous-growth model
//' sim$border_growth_model <- FALSE
//'
//' # now it is set to `FALSE`
//' sim$border_growth_model
//'
//' # switch back to the border-growth model
//' sim$border_growth_model <- FALSE
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .property("border_growth_model", &TissueSimulation::is_border_growth_model,
                  &TissueSimulation::set_border_growth_model,
                  "Switch between homogeneous and border driven growth models.")

//' @name TissueSimulation$get_clock
//' @title Getting the simulated time
//' @description This method returns the current simulation time.
//' @return The time simulated by the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to time 40
//' sim$run_up_to_time(40)
//'
//' # get the simulated time
//' sim$get_clock()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_clock", &TissueSimulation::get_clock,
                "Get the current simulation time")

//' @name TissueSimulation$get_cell
//' @title Getting one of the tissue cells
//' @description This method collects some data of the aimed cell without
//'   altering the tissue.
//' @param x The position of the aimed cell on the x axis.
//' @param y The position of the aimed cell on the y axis.
//' @return A data frame reporting `cell_id`, `mutant`, `position_x`, and
//'   `position_y` of the aimed cell. If the simulation has epigenetic
//'   states, the data frame also contains the column `epistate`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # collect a cell in the tissue
//' sim$get_cell(501, 502)
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_cell",
                (List (TissueSimulation::*)(const RE::AxisPosition &,
                                            const RE::AxisPosition &)
                     const)(&TissueSimulation::get_cell),
                "Get one cell from the simulated tissue")

//' @name TissueSimulation$get_cells
//' @title Getting the tissue cells
//' @description This method returns information about tumour tissue cells
//' @details It collects some data about the cells in the tissue
//'   without altering the tissue itself. The pairs of optional parameters
//'   `lower_corner` and `upper_corner` define a frame of the tissue in
//'   which the data are sampled. The optional parameters `mutant_filter`
//'   and `epigenetic_filter` filter the collected cell data according to
//'   the cell mutant and epigenetic state.
//' @param lower_corner The lower-left corner of the selection frame
//'   (optional).
//' @param upper_corner The upper-right corner of the selection frame
//'   (optional).
//' @param mutant_filter The vector of the to-be-selected mutant names
//'   (optional).
//' @param epigenetic_filter The vector of the to-be-selected epigenetic states
//'   (optional).
//' @return A data frame having the columns `cell_id`, `mutant`,
//'   `position_x`, `position_y`, and, when the simulation has epigenetic
//'   states, `epistate`. Each row in the data frame describes a cell
//'   that satisfies the provided filters and lays in the specified frame.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' sim$death_activation_level <- 100
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.3, death = 0.1, E2 = 0.02),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 20)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # collect all cells in the tissue
//' cells <- sim$get_cells()
//'
//' # print some of them
//' head(cells)
//'
//' # let us define a function to print some statistics
//' print_statistics <- function(cells) {
//'   cat(paste("num of cells:", nrow(cells)))
//'   cat(paste("\nmutants:", paste(unique(cells$mutant), collapse = " ")))
//'   cat(paste("\nepigenetic states:", paste(unique(cells$epistate),
//'                                           collapse = " ")))
//'   cat(paste0("\nframe: [",
//'              paste(range(cells$position_x), collapse = ","), "]x[",
//'              paste(range(cells$position_y), collapse = ","), "]\n\n"))
//' }
//'
//' print_statistics(cells)
//'
//' # get the cells in the frame [495,505]x[490,500]
//' cells <- sim$get_cells(lower_corner = c(495, 495),
//'                        upper_corner = c(505, 505))
//'
//' print_statistics(cells)
//'
//' # cells can be filtered by mutant name...
//' cells <- sim$get_cells(mutant_filter = c("A"),
//'                        epigenetic_filter = c("E1", "E2"))
//'
//' print_statistics(cells)
//'
//' # ..., by epigenetic state, ...
//' cells <- sim$get_cells(mutant_filter = c("A", "B"),
//'                        epigenetic_filter = c("E2"))
//'
//' print_statistics(cells)
//'
//' # ..., or by position, mutant, and epigenetic state.
//' cells <- sim$get_cells(lower_corner = c(495, 495), upper_corner = c(505, 505),
//'                        mutant_filter = c("A", "B"),
//'                        epigenetic_filter = c("E2"))
//'
//' print_statistics(cells)
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_cells",
                (List (TissueSimulation::*)(const std::vector<RE::AxisPosition> &,
                                            const std::vector<RE::AxisPosition> &,
                                            const std::vector<std::string> &,
                                            const std::vector<std::string> &)
                     const)(&TissueSimulation::get_cells),
                "Get cells from the simulated tissue")
        .method("get_cells",
                (List (TissueSimulation::*)(const SEXP &, const SEXP &)
                     const)(&TissueSimulation::get_cells),
                "Get cells from the simulated tissue")
        .method("get_cells",
                (List (TissueSimulation::*)(const SEXP &)
                     const)(&TissueSimulation::get_cells),
                "Get cells from what was the simulated tissue status before a sampling")
        .method("get_cells",
                (List (TissueSimulation::*)() const)(&TissueSimulation::get_cells),
                "Get cells from the simulated tissue")

//' @name TissueSimulation$get_name
//' @title Getting the simulation name
//' @description This method returns the simulation name
//' @return The simulation name, which corresponds to the name of the directory
//'   in which the simulation is saving its progresses.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # Expecting "test"
//' sim$get_name()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_name", &TissueSimulation::get_name, "Get the simulation name")

//' @name TissueSimulation$get_lineage_graph
//' @title Getting the simulation lineage graph
//' @description This method returns the lineage graph of the simulation.
//' @details At the beginning of the computation only the species of the
//'   added cells are present in the tissue. As the simulation proceeds new
//'   species arise as a consequence of either mutant mutations or epigenetic
//'   switches. The *lineage graph* stores these species evolutions and it
//'   reports the first occurrence time of any species-to-species transition.
//' @return A data frame reporting `ancestor`, `progeny`, and `first_occurrence`
//'   of each species-to-species transition.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.3, death = 0.1, E2 = 0.02),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 20)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' sim$get_lineage_graph()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_lineage_graph", &TissueSimulation::get_lineage_graph,
                "Get the simulation lineage graph")

//' @name TissueSimulation$get_tissue_size
//' @title Getting the simulated tissue size
//' @description This method returns the size of the simulated tissue.
//' @return The vector `c(x_size, y_size)` of the simulated tissue.
//' @examples
//' # create a simulation having size 1200x900
//' sim <- TissueSimulation(width = 1200, height = 900)
//'
//' # get the tissue size, i.e., expecting c(1200,900)
//' sim$get_tissue_size()
//' @seealso [TissueSimulation()], <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_tissue_size", &TissueSimulation::get_tissue_size,
                "Get the simulation tissue size")

//' @name TissueSimulation$get_added_cells
//' @title Getting the cells manually added to the simulation
//' @description This method returns the cells manually added to
//'   the simulation.
//' @return A data frame having the columns `mutant`, `position_x`,
//'   `position_y`, `time`, and, when the simulation has epigenetic
//'   states, `epistate`. The data frame contains a row for each
//'   cell manually added to the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # get the cells
//' sim$get_added_cells()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_added_cells", &TissueSimulation::get_added_cells,
                "Get the cells manually added to the simulation")

//' @name TissueSimulation$get_counts
//' @title Counting the cell number
//' @description This method returns the current number of cells per
//'   species and that since the simulation began.
//' @return A data frame reporting `mutant`, `counts`, `overall`, and,
//'   when the simulation has epigenetic states, `epistate`. The data
//'   frame contains a row for each species in the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", c(duplication = 0.4, death = 0.01))
//'
//' # schedule an evolution from mutant "A" to mutant "B" at time 20
//' sim$schedule_mutation("A", "B", 20)
//'
//' # place a cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.4, death = 0.1, E2 = 0.02),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 10)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_counts", &TissueSimulation::get_counts,
                "Get the current number of cells and that "
                "throughout the entire simulation")

//' @name TissueSimulation$get_count_history
//' @title Getting the history of the number of cells per species
//' @description This method returns a data frame reporting the number of
//'   species cells in each sampled simulation time.
//' @return A data frame reporting `mutant`, `counts`, and `time` for each
//'   species, and for each sampled time. When the simulation has epigenetic
//'   states, the data frame also contains the column `epistate`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", c(duplication = 0.15, death = 0.05))
//'
//' # schedule an evolution from mutant "A" to mutant "B" at time 10
//' sim$schedule_mutation("A", "B", 10)
//'
//' # place a cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # set delta time between species counting to 20
//' sim$history_delta <- 20
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # get the history of species counts
//' sim$get_count_history()
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.4, death = 0.1, E2 = 0.02),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 10)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # set delta time between species counting to 30
//' sim$history_delta <- 30
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_count_history()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_count_history",
                (List (TissueSimulation::*)() const) &
                    TissueSimulation::get_count_history,
                "Get the number of simulated events per species along the computation")

//' @name TissueSimulation$get_firings
//' @title Getting the number of fired events
//' @description This method returns a data frame reporting the current
//'   number of simulated events per species.
//' @return A data frame having the `event`, `mutant`, `fired`, and, when
//'   the simulation has epigenetic states, `epistate`. Each row reports
//'   event of a given type have been fired in a species.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", c(duplication = 0.15, death = 0.05))
//'
//' # schedule an evolution from mutant "A" to mutant "B" at time 10
//' sim$schedule_mutation("A", "B", 10)
//'
//' # place a cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # get the history of species counts
//' sim$get_firings()
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.4, death = 0.1, E2 = 0.02),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # schedule a mutation from "A" to "B"
//' sim$schedule_mutation("A", "B", 10)
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_firings()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_firings", &TissueSimulation::get_firings,
                "Get the current number of simulated events per species")

//' @name TissueSimulation$get_firing_history
//' @title Getting the fired event history
//' @description This method returns a data frame reporting the number of
//'   events fired up to each sampled simulation time.
//' @return A data frame reporting `event`, `mutant`, `fired`, and `time`
//'   for each event type, for each species, and for each sampled time.
//'   Whenever, the simulation has epigenetic states, the data frame also
//'   contains the column `epistate`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # set delta time between species counting to 30
//' sim$history_delta <- 30
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # get the number of event fired per event and species
//' sim$get_firing_history()
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_firing_history",
                (List (TissueSimulation::*)() const) &
                    TissueSimulation::get_firing_history,
                "Get the number of simulated events per species along the computation")

//' @name TissueSimulation$get_rates
//' @title Getting the simulation rates
//' @description This method returns the rates of the simulation.
//' @details This method returns a data frame containing the simulation rates.
//'    A rate is not included in the returned data frame if and only if it was not
//'    set during the system specification. In these cases, the rate is assumed to
//'    be 0 by default.
//' @param complete A Boolean flag to get also the rates that have not been set
//'    (default: `FALSE`).
//' @return A data frame containing the simulation rates. If the
//'    simulation has epigenetic states, the data frame has 5 columns: `mutant`,
//'    `epistate`, `event`, `first.child.epistate`, and `rate`. The columns `mutant`
//'    and `epistate` store the mutant and the epigenetic state of the cell from
//'    which the event may occur. The columns `event` and `rate` maintain the name
//'    and the rate of the event. Finally, the column `first.child.epistate` reports
//'    the epigenetic state of potential first child due to the event. For instance,
//'    when the event is `duplication`, the first child has the same epigenetic state
//'    of its parent. Instead, when the event is `switch`, the column
//'    `first.child.epistate` contains an epigenetic state diffent from that of the
//'    origin cell. In case of the event `death`, the column `first.child.epistate`
//'    is set to NA.
//'    When the simulation has no epigenetic states, the returned data frame
//'    exclusively contains the columns `mutant`, `event`, and `rate`.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add a mutant
//' sim$add_mutant("A")
//'
//' # set its duplication rate
//' sim$set_rate("A", "duplication", 0.1)
//'
//' # get the rates that have been set
//' sim$get_rates()
//'
//' # get all simulation rates
//' sim$get_rates(TRUE)
//'
//' # add epigenetic states
//' sim$add_epigenetic_states(c("E1", "E2"))
//'
//' # set some of the rates of "A[E1]" and "A[E2]"
//' sim$set_rate("A[E1]", "duplication", 0.1)
//' sim$set_rate("A[E2]", "duplication", 0.1)
//' sim$set_rate("A[E1]", "death", 0.2)
//'
//' # get the rates that have been set
//' sim$get_rates()
//'
//' # get all simulation rates
//' sim$get_rates(TRUE)
//' @seealso [TissueSimulation$set_rate()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_rates",
                (List (TissueSimulation::*)() const)(
                    &TissueSimulation::get_rates),
                "Get the rates of the simulation")
        .method("get_rates",
                (List (TissueSimulation::*)(const SEXP&) const)(
                    &TissueSimulation::get_rates),
                "Get the rates of the simulation")

//' @name TissueSimulation$get_rates_update_history
//' @title Retrieving the rates update history
//' @description This method retrieves the simulation rates
//'   update history.
//' @return A data frame containing the set simulation rates. If the simulation has
//'    epigenetic states, the data frame has 6 columns: `time`, `mutant`, `epistate`,
//'    `event`, `first.child.epistate`, and `rate`. The column `time` contains the
//'    rate setting time. The columns `mutant` and `epistate` store the mutant and
//'    the epigenetic state of the cell from which the event may occur.
//'    The columns `event` and `rate` maintain the name and the rate of the event.
//'    Finally, the column `first.child.epistate` reports the epigenetic state of
//'    potential first child due to the event. For instance, when the event is
//'    `duplication`, the first child has the same epigenetic state of its parent.
//'    Instead, when the event is `switch`, the column `first.child.epistate`
//'    contains an epigenetic state diffent from that of the origin cell. In case
//'    of the event `death`, the column `first.child.epistate` is set to NA.
//'    When the simulation has no epigenetic states, the returned data frame
//'    exclusively contains the columns `time`, `mutant`, `event`, and `rate`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # place a cell of "A"
//' sim$place_cell("A", 500, 500)
//'
//' # set the death activation level
//' sim$death_activation_level <- 100
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # set the death rate of "A" to 0.9
//' sim$set_rate("A", "death", 0.9)
//'
//' # we changed our mind *before* running the simulation
//' # and we reset the death rate to 0.05
//' sim$set_rate("A", "death", 0.05)
//'
//' # simulate up to time 80
//' sim$run_up_to_time(80)
//'
//' # set the death rate to 0.5
//' sim$set_rate("A", "death", 0.5)
//'
//' # simulate up to time 80+1
//' sim$run_up_to_time(sim$get_clock()+1)
//'
//' # get the rates update history
//' sim$get_rates_update_history()
//' @seealso [TissueSimulation$set_rates()],
//'   [TissueSimulation$get_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_rates_update_history", &TissueSimulation::get_rates_update_history,
                "Get the rates update history")

//' @name TissueSimulation$get_samples_info
//' @title Retrieving sample information
//' @description This method retrieves information about
//'   the samples collected along the simulation.
//' @return A data frame containing, for each sample collected
//'   during the simulation, the columns `name`, `time`, `id`,
//'   `ymin`, `xmin`, `ymax`, `xmax`, `tumour_cells`, and
//'   `tumour_cells_in_bbox`. The columns `ymin`, `xmin`, `ymax`,
//'   `xmax` report the boundaries of the sample bounding box, while
//'   `tumour_cells` and `tumour_cells_in_bbox` are the number of tumour
//'   cells in the sample and in the bounding box, respectively.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", list(duplication = 0.3, death = 0.01))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue until "A" consists of 50000 cells
//' sim$run_up_to_size("A", 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner = c(450, 475),
//'                  upper_corner = c(500, 550))
//'
//' # simulate 1 time unit more
//' sim$run_up_to_time(sim$get_clock()+1)
//'
//' # sample the region [500,520]x[525,550]
//' sim$sample_cells("S2", lower_corner = c(500, 525),
//'                  upper_corner = c(520, 550))
//'
//' # get information about all the collected
//' # samples, i.e, S1 and S2
//' sim$get_samples_info()
//' @seealso [TissueSimulation$sample_cells()],
//'   [SampleForest$get_samples_info()],
//'   [PhylogeneticForest$get_samples_info()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("get_samples_info",
                (List (TissueSimulation::*)() const) & TissueSimulation::get_samples_info,
                "Get some pieces of information about the collected samples")

//' @name TissueSimulation$history_delta
//' @title The delta time between time series samples
//' @description The maximum time between two consecutive
//'   time series data samples.
//' @details This property is the maximum time between two consecutive
//'   time series data samples. Notice that this property differs from
//'   <code>[TissueSimulation$snapshot_triggers](TissueSimulation-cash-snapshot_triggers.md)</code>.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # get the delta time between two time series samples (0 by default)
//' sim$history_delta
//'
//' # set the delta time between two time series samples
//' sim$history_delta <- 20
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>,
//'   <code>[TissueSimulation$snapshot_triggers](TissueSimulation-cash-snapshot_triggers.md)</code>
        .property("history_delta", &TissueSimulation::get_history_delta,
                  &TissueSimulation::set_history_delta,
                  "The sampling delta for the get_*_history functions")

//' @name TissueSimulation$mutate_progeny
//' @title Generating a mutated progeny
//' @description This method generates a mutated progeny.
//' @details It simulates both the duplication of the cell in the
//'   specified position and the birth of one cells of a given
//'   mutant that preserves the epigenetic status of the original cell.
//'   The mutated cell will be located in the position of its parent.
//' @param cell_position The position of the cell whose offspring will mutate.
//' @param mutated_mutant The mutant of the mutated cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # add mutant "B" and set its species rates
//' sim$add_mutant("B",
//'                list(E1 = list(duplication = 0.15, death = 0.3, E2 = 0.1),
//'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation up to time 70
//' sim$run_up_to_time(70)
//'
//' # get the number of cells per species. No cell in "B" yet.
//' sim$get_counts()
//'
//' # duplicate the cell in position (503, 492). One of
//' # its direct descendants will have mutant "B"
//' # sim$mutate_progeny(503, 492, "B")
//'
//' # the output of `choose_cell_in`, `choose_border_cell_in` and `get_cell`
//' # can also be used as input for `mutate_progeny`
//' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
//'
//' # get the number of cells per species again.
//' # Now, "B" consists of one cell
//' sim$get_counts()
//' @seealso [TissueSimulation$choose_cell_in()],
//'   [TissueSimulation$choose_border_cell_in()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("mutate_progeny",
                (void (TissueSimulation::*)(const List &, const std::string &))(
                    &TissueSimulation::mutate_progeny),
                "Duplicate a cell and mutate one of its children")
        .method("mutate_progeny",
                (void (TissueSimulation::*)(
                    const RE::AxisPosition &, const RE::AxisPosition &,
                    const std::string &))(&TissueSimulation::mutate_progeny),
                "Duplicate a cell and mutate one of its children")

//' @name TissueSimulation$run_up_to_time
//' @title Simulating cell evolution
//' @param time The final simulation time.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation without epigenetic states
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", list(duplication = 0.3, death = 0.01))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue up to simulate timed 40
//' sim$run_up_to_time(40)
//' @seealso [TissueSimulation$run_up_to_event()],
//'    [TissueSimulation$run_up_to_size()], [TissueSimulation$run_until()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("run_up_to_time",
                (void (TissueSimulation::*)(const CLONES::Time &,
                                            const bool))&TissueSimulation::run_up_to_time,
                "Simulating the system up to the specified simulation time")
        .method("run_up_to_time",
                (void (TissueSimulation::*)(
                    const CLONES::Time &))&TissueSimulation::run_up_to_time,
                "Simulating the system up to the specified simulation time")

//' @name TissueSimulation$run_up_to_event
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until the number of events
//'   that have occurred to cells of a species reaches a specified threshold.
//' @param event The considered event, i.e., `growth`, `death`, or `switch`.
//' @param species The species whose event number is considered.
//' @param num_of_events The threshold for the event number.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # simulate the cell evolution until the number of epigenetic events from
//' # the species "A[E2]" is less than 100.
//' sim$run_up_to_event(event = "switch", species = "A[E2]",
//'                     num_of_events = 100)
//'
//' sim
//' @seealso [TissueSimulation$run_up_to_time()],
//'    [TissueSimulation$run_up_to_size()], [TissueSimulation$run_until()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("run_up_to_event",
                (void (TissueSimulation::*)(
                    const std::string &, const std::string &, const size_t &,
                    const bool))&TissueSimulation::run_up_to_event,
                "Simulating the system up to the specified number of events")
        .method("run_up_to_event",
                (void (TissueSimulation::*)(
                    const std::string &, const std::string &,
                    const size_t &))&TissueSimulation::run_up_to_event,
                "Simulating the system up to the specified number of events")

//' @name TissueSimulation$run_up_to_size
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until the number of
//'   cells in a species reaches a specified threshold.
//' @param species The species whose number of cells is considered.
//' @param num_of_cells The threshold for the cell number.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # simulate the tissue until the species "A[E2]" account for 100
//' # contemporary cells
//' sim$run_up_to_size(species = "A[E2]", num_of_cells = 100)
//'
//' sim$get_counts()
//' @seealso [TissueSimulation$run_up_to_time()],
//'    [TissueSimulation$run_up_to_event()], [TissueSimulation$run_until()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "run_up_to_size",
            (void (TissueSimulation::*)(const std::string &, const size_t &,
                                        const bool))&TissueSimulation::run_up_to_size,
            "Simulating the system up to the specified number of cells in the species")
        .method(
            "run_up_to_size",
            (void (TissueSimulation::*)(const std::string &,
                                        const size_t &))&TissueSimulation::run_up_to_size,
            "Simulating the system up to the specified number of cells in the species")

//' @name TissueSimulation$run_until
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until a formula does not
//'    hold.
//' @param formula The formula that will be satisfied at the end of the
//'    simulation.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place an "A[E1]" cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # get the variable representing the simulation time
//' v_time <- sim$var("Time")
//'
//' # get the variable representing the cardinality of A[E1]
//' va_e1 <- sim$var("A[E1]")
//'
//' # get the variable representing the cardinality of A[E2]
//' va_e2 <- sim$var("A[E2]")
//'
//' # get the variable representing the number of epigenetic
//' # switches from A[E1]
//' va_ps <- sim$var("A[E1].switches")
//'
//' # build a condition stating that the cardinality of A[E1]
//' # doubles that of A[E2]
//' c1 <- va_e1 >= 2*va_e2
//'
//' # build a condition that holds when there are more than
//' # 100000 live cells of mutant A
//' c2 <- va_e1 + va_e2 > 1e5
//'
//' # build a condition that holds when less than 4000 switched
//' # from A[E1] have occurred
//' c3 <- va_ps < 4000
//'
//' # build a condition that holds when 40 time unit have been
//' # simulated at least
//' c4 <- v_time >= 40
//'
//' # build a condition that holds when c4 and at least one
//' # among c1, c2, and c3 hold
//' c5 <- c4 & (c1 | c2 | c3)
//' c5
//'
//' # run the simulation while c5 does not hold
//' sim$run_until(c5)
//'
//' sim
//' @seealso [TissueSimulation$var()],
//'   [TissueSimulation$run_up_to_time()], [TissueSimulation$run_up_to_event()],
//'   [TissueSimulation$run_up_to_size()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("run_until",
                (void (TissueSimulation::*)(const Logics::Formula &,
                                            const bool))&TissueSimulation::run_until,
                "Simulating the system until a formula is *not* satisfied")
        .method("run_until",
                (void (TissueSimulation::*)(
                    const Logics::Formula &))&TissueSimulation::run_until,
                "Simulating the system until a formula is *not* satisfied")

//' @name TissueSimulation$sample_cells
//' @title Sampling a set of cells
//' @description This method samples a set of tumour cells.
//' @details It removes the cells from the simulated tissue and
//'   stores them in a sample that can be subsequently
//'   retrieved to build a sample forest.
//' @param sample_name The name of the sample.
//' @param lower_corner The lower corner of the sample bounding box (optional
//'   in pair with `upper_corner`).
//' @param upper_corner The upper corner of the sample bounding box (optional
//'   in pair with `lower_corner`).
//' @param num_of_cells The maximum number of tumour cells to collect
//'   (optional).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", list(duplication = 0.2, death = 0.01))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue until "A" consists of 50000 cells
//' sim$run_up_to_size("A", 50000)
//'
//' # randomly sample 50 tumour cells from the tissue
//' sim$sample_cells(sample_name = "S1", num_of_cells = 50)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells(sample_name = "S2",
//'                  lower_corner = c(450, 475), upper_corner = c(500, 550))
//'
//' # randomly sample 50 tumour cells from the tissue region [500,550]x[500,550]
//' sim$sample_cells(sample_name = "S3",
//'                  lower_corner = c(500, 500), upper_corner = c(550, 550),
//'                  num_of_cells = 50)
//'
//' sim$get_samples_info()
//' @seealso [TissueSimulation$get_samples_info()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method(
            "sample_cells",
            (void (TissueSimulation::*)(
                const std::string &,
                const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
                const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner,
                const size_t &num_of_cells) const)(&TissueSimulation::sample_cells),
            "Sample a rectangular region of the tissue")
        .method(
            "sample_cells",
            (void (TissueSimulation::*)(
                const std::string &,
                const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
                const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner)
                 const)(&TissueSimulation::sample_cells),
            "Sample a rectangular region of the tissue")
        .method(
            "sample_cells",
            (void (TissueSimulation::*)(const std::string &, const size_t &num_of_cells)
                 const)(&TissueSimulation::sample_cells),
            "Sample a rectangular region of the tissue")

//' @name TissueSimulation$set_rate
//' @title Set the rate of an event
//' @description This method sets the species' rate of an event.
//' @param species The species whose rate must be set.
//' @param event_name The name of the event whose rate should be set.
//' @param dest Either the species or the epigenetic state of
//'   one of the children due to the event (to be specified for the
//'   epigenetic switch event only).
//' @param rate The rate of the event.
//' @examples
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2", "E3"))
//'
//' # add mutant "A"
//' sim$add_mutant("A")
//'
//' # set the duplication and death rates of the species "A[E1]"
//' sim$set_rate("A[E1]", "duplication", 0.1)
//' sim$set_rate("A[E1]", "death", 0.2)
//'
//' # setting the switch event rates from "A[E1]" to "A[E2]" and "A[E3]"
//' sim$set_rate("A[E1]", "switch", "E2", 0.0001)
//' sim$set_rate("A[E1]", "switch", "A[E3]", 0.0002)
//'
//' sim$get_rates()
//' @seealso [TissueSimulation$set_rates()], [TissueSimulation$get_rates()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("set_rate", (void (TissueSimulation::*)(const std::string &,
                    const std::string &, const double &))(&TissueSimulation::set_rate),
                "Set the rate of an event")
        .method("set_rate", (void (TissueSimulation::*)(const std::string &,
                    const std::string &, std::string,
                    const double &))(&TissueSimulation::set_rate),
                "Set the rate of an event")

//' @name TissueSimulation$set_rates
//' @title Set the tissue simulation rates
//' @description This method sets tissue simulation rates.
//' @details This method can set the rates of multiple species. It accepts as
//'   parameters both a list of a data frame specifying the rates to be set.
//'
//'    - The list must be a named list whose names correspond to the species
//'      whose rates must be set. Each element in the list represents the rates
//'      to be set for the corresponding species and is itself a named list
//'      whose names must belong to the set containing `"death"`,
//'      `"duplication"`, and, if the simulation includes epigenetic states,
//'      the names of the known epigenetic states. The values of the element
//'      named `"duplication”` and `“death”` represent the new species’
//'      duplication and death rates, respectively. Instead, the values whose
//'      names are among the known epigenetic states indicate the new switch
//'      rate to the specified epigenetic state (see below for examples).
//'    - The data frame must contain at least three columns: `mutant`, `event`,
//'      and `rate`. If the simulation includes epigenetic states, the data
//'      frame must also include the columns `epistate` and
//'      `first.child.epistate`. Each row in the data frame declares the new
//'      rate value for an event in a species. The columns `mutant` and `event`
//'      represent the species mutants and event types. The values in the
//'      former column must be known mutant names, while those in the latter
//'      must be among `“death”`, `“duplication”`, and, when the simulation
//'      includes epigenetic states, `“switch”`. The columns
//'      `epistate` and `first.child.epistate` denote the epigenetic state of
//'      the species whose rate is to be set and the epigenetic state of the
//'      first child due to the event. When the event is `“duplication”`, the
//'      `first.child.epistate` value should be equal to the value in the
//'      column `epistate`; when the event is `“switch”`, it
//'      must be a known epigenetic state different from that reported in the
//'      column epistate; when the event is `“death”`, it can be `NA`. Finally,
//'      the column `rate` contains the new rate values and must be numerical.
//'
//'   Notice that no new epigenetic state nor mutants will be added to the
//'   simulation by this method. Any mention to a non-existant mutant or
//'   epigenetic state ends the execution with an error.
//'
//' @param rates Either a list or a data frame of the rates to be set.
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutants "A", "B", and "C"
//' sim$add_mutants(c("A", "B", "C"))
//'
//' # no rate is reported because none has been set yet
//' sim$get_rates()
//'
//' # setting "A"'s duplication and death rates and "B"'s death rate
//' # All remaining rates are set to zero by default
//' sim$set_rates(list(A = list(duplication = 0.3, death = 0.2),
//'                    B = list(death = 0.2)))
//'
//' # get_rates() only reports the rates that have been explicitly set
//' sim$get_rates()
//'
//' # adding epigenetic states
//' sim$add_epigenetic_states(c("E1", "E2", "E3"))
//'
//' # species before adding epigenetic states do not exists anymore
//' # and no rates have been set for the new species
//' sim$get_rates()
//'
//' # set some rates for the new species
//' sim$set_rates(list("A[E1]" = list(duplication = 0.3, death = 0.2,
//'                                   E2 = 0.01, "A[E3]" = 0.04),
//'                    "B[E2]" = list(death = 0.2, E3 = 0.01)))
//'
//' sim$get_rates()
//' @examples
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutants "A", "B", and "C"
//' sim$add_mutants(c("A", "B", "C"))
//'
//' # build a data frame for setting "A"'s duplication and death
//' # rates and "B"'s death rate as well
//' df_rates <- data.frame(
//'  mutant = c("A", "A", "B"),
//'  event = c("duplication", "death", "duplication"),
//'  rate = c(0.3, 0.2, 0.2)
//' )
//'
//' df_rates
//'
//' # setting the rates by using the data frame
//' sim$set_rates(df_rates)
//'
//' # get_rates() only reports the set rates.
//' sim$get_rates()
//'
//' # adding epigenetic states
//' sim$add_epigenetic_states(c("E1", "E2", "E3"))
//'
//' # now we need to add the "epistate" and `first.child.epistate`
//' # columns to the data frame
//' df_rates[["epistate"]] <- c("E1", "E1", "E2")
//' df_rates[["first.child.epistate"]] <- c("E1", NA, "E2")
//'
//' # load dplyr to simplify the next part of the example
//' library(dplyr)
//'
//' # we also may set some switch rates
//' df_rates <- df_rates %>%
//'   add_row(mutant = "A", epistate = "E1", event = "switch",
//'           first.child.epistate = "E2", rate = 0.01) %>%
//'   add_row(mutant = "A", epistate = "E1", event = "switch",
//'           first.child.epistate = "E3", rate = 0.04) %>%
//'   add_row(mutant = "B", epistate = "E2", event = "switch",
//'           first.child.epistate = "E3", rate = 0.01)
//'
//' df_rates
//'
//' sim$set_rates(df_rates)
//'
//' # get_rates() only reports the set rates. All remaining rates are
//' # set to zero by default
//' sim$get_rates()
//' @seealso [TissueSimulation$set_rate()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("set_rates", (void (TissueSimulation::*)(
                    const Rcpp::List& rate_list))(&TissueSimulation::set_rates),
                "Set tissue simulation rates")

//' @name TissueSimulation$search_sample
//' @title Searching for a rectangular tissue sample
//' @description This method searches a rectangular tissue sample.
//' @details The aimed sample must satisfy the specified number of cells.
//'   The size of the sample is also provided a parameter of the method.
//'   The complexity of this method is
//'   \eqn{O(|\textrm{tissue width}|*|\textrm{tissue height}|)}.
//' @param min_num_of_cells A named integer vector reporting the minimum number
//'   of cells per species or mutant.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @return If a rectangular sample satisfying the provided constraints can
//'   be found, the corresponding rectangle.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", list(duplication = 0.2, death = 0.01))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue until "A" consists of 100 cells
//' sim$run_up_to_size("A", 100)
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", list(duplication = 0.3, death = 0.01))
//'
//' # mutate a border cell in "A" into "B"
//' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
//'
//' # run the simulation until "B" consists of 1000 cells
//' sim$run_up_to_size("B", 1000)
//'
//' # find a 50x50 sample containing 80 "B" cells and 10 "A" cells at least
//' sim$search_sample(c("A" = 10, "B" = 80), 50, 50)
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("search_sample", &TissueSimulation::search_sample,
                "Search a rectangular sample containing a given number of cells")

//' @name TissueSimulation$search_samples
//' @title Searching rectangular tissue samples
//' @description This method searches a set of rectangular tissue samples.
//' @details The aimed samples mush satisfy the specified number of cells.
//'   The sizes of the samples are also provided a parameter of the method.
//'   This method takes asymptotic time
//'   \eqn{O(|\textrm{tissue width}|*|\textrm{tissue height}|)}.
//' @param min_num_of_cells A named integer vector reporting the minimum number
//'   of cells per species or mutant.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @param n_samples The number of searched samples.
//' @param seed The seed of the random generator the select the samples
//'   among those satisfying the constraints (optional).
//' @return A vector of `n_samples` rectangular tissue samples that
//'   satisfy the aimed constraints.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation(width = 150, height = 150)
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", list(duplication = 0.20, death = 0.01))
//'
//' # place an "A" cell in the tissue
//' sim$place_cell("A", 75, 75)
//'
//' # simulate the tissue until "A" consists of 300 cells
//' sim$run_up_to_size("A", 130)
//'
//' # add mutant "B" and set its rates
//' sim$add_mutant("B", list(duplication = 0.3, death = 0.01))
//'
//' # mutate a border cell in "A" into "B"
//' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
//'
//' # run the simulation until "B" consists of 1000 cells
//' sim$run_up_to_size("B", 4000)
//'
//' # plot the tissue as it is
//' plot <- plot_tissue(sim)
//'
//' # find 3 50x50 samples containing 80 "B" cells and 100 "A" cells
//' # at least
//' bboxes <- sim$search_samples(min_num_of_cells = c("A" = 100, "B" = 80),
//'                              width = 25, height = 25,
//'                              n_samples = 3)
//' bboxes
//'
//' # plot the found bounding boxes
//' for (bbox in bboxes) {
//'   plot <- plot +
//'     ggplot2::geom_rect(xmin = bbox$lower_corner[1],
//'                        xmax = bbox$upper_corner[1],
//'                        ymin = bbox$lower_corner[2],
//'                        ymax = bbox$upper_corner[2],
//'                        fill = NA, color = "black")
//' }
//'
//' plot
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("search_samples",
                (std::vector<TissueRectangle> (TissueSimulation::*)(
                    const Rcpp::IntegerVector &, const uint16_t &, const uint16_t &,
                    const size_t) const)(&TissueSimulation::search_samples),
                "Search rectangular samples containing a given number of cells")
        .method("search_samples",
                (std::vector<TissueRectangle> (TissueSimulation::*)(
                    const Rcpp::IntegerVector &, const uint16_t &, const uint16_t &,
                    const size_t, int) const)(&TissueSimulation::search_samples),
                "Search rectangular samples containing a given number of cells")

//' @name TissueSimulation$var
//' @title Building a simulation status variable
//' @description This method builds a simulation status variable.
//' @details This method builds a logic variable representing one of the
//'   simulation quantities among:
//'   -  cardinality of a species
//'   -  number of event among duplications, deaths, and epigenetic switches
//'   -  elapsed evolution time
//' @param variable_description The description of the variable to be built.
//'   When `variable_description` is the string `"Time"`, the elapsed
//'   simulation time variable is returned. If `variable_description` is
//'   set to a species name, then the variable representing the cardinality
//'   of the species is built. Finally, when the parameter is a species name
//'   followed by `.` and one among `duplications`, `deaths`, or `switches`,
//'   the variable representing the number of event of the specified type
//'   occurred since the computation beginning in the species.
//' @return A variable representing the simulation quantity according to
//'   the parameter `variable_description`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.01, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # get the variable representing the simulation time
//' sim$var("Time")
//'
//' # get the variable representing the cardinality of A[E1]
//' sim$var("A[E1]")
//'
//' # get the variable representing the cardinality of A[E2]
//' sim$var("A[E2]")
//'
//' # get the variable representing the number of duplications
//' # from A[E2]
//' sim$var("A[E2].duplications")
//'
//' # get the variable representing the number of epigenetic
//' # switches from A[E1]
//' sim$var("A[E1].switches")
//' @seealso [TissueSimulation$run_until()],
//'   <code>[TissueSimulation](TissueSimulation_class.md)</code>
        .method("var", &TissueSimulation::get_var,
                "Get a variable representing a simulation quantity");

//' @name recover_simulation
//' @title Loading a simulation
//' @description This method loads a simulation from the disk.
//' @usage recover_simulation(path)
//' @param path The path of the directory in which the simulation
//'   was saved or that of a snapshot file.
//' @examples
//' # set the random seed for repeatability
//' set.seed(0)
//'
//' # create a simulation having name "recover_simulation_test" and
//' # save its snapshots in a local directory
//' sim <- TissueSimulation("recover_simulation_test",
//'                         epigenetic_states = c("E1", "E2"),
//'                         save_directory = TRUE)
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
//'                     E2 = list(duplication = 0.01, death = 0.005, E1 = 0.01)))
//'
//' # place a cell in the tissue
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # simulate up to time 50
//' sim$run_up_to_time(50)
//'
//' # simulate up to time 75. Each new simulation step produces
//' # a new snapshot
//' sim$run_up_to_time(75)
//'
//' # simulate up to time 80
//' sim$run_up_to_time(80)
//'
//' # show the simulation
//' sim
//'
//' # remove the object sim from the environment
//' rm(list = c("sim"))
//'
//' # the object pointed by sim does not exist any more
//' exists("sim")
//'
//' # recover the simulation from the directory "recover_simulation_test"
//' sim <- recover_simulation("recover_simulation_test")
//'
//' sim
//'
//' # get the sim clock
//' sim$get_clock()
//'
//' # get the snapshot information
//' snapshot_info <- sim$get_snapshot_info()
//'
//' # load the second snapshot
//' sim <- recover_simulation(snapshot_info[["file"]][2])
//'
//' sim
//'
//' # get the sim clock
//' sim$get_clock()
//'
//' # delete dump directory
//' unlink("recover_simulation_test", recursive = TRUE)
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
    function("recover_simulation", &TissueSimulation::load, "Recover a simulation");

//' @name TissueSimulation
//' @title Building a new simulation
//' @description This method builds a new simulation.
//' @usage TissueSimulation(name, width = 1000, height = 1000, save_directory = FALSE)
//' @param name The name of the simulation (default:
//'   "`ProCESS_<year>_<hour><minute><second>`").
//' @param width The width of the simulated tissue (default: 1000).
//' @param height The height of the simulated tissue (default: 1000).
//' @param save_directory A flag to save simulation data in the working directory
//'   (default: `FALSE`).
//' @param rates A data frame specifying the simulation species and their rates
//'   (default: `NULL`). See [TissueSimulation$set_rates()] for the data frame
//'   specification. Differently from [TissueSimulation$set_rates()] the species
//'   are automatically added to the new simulation.
//' @param epigenetic_states A list of epigenetic states
//'   (default: `NULL`).
//' @param seed The seed for the pseudo-random generator (optional).
//' @examples
//' # create a TissueSimulation object storing binary dump in a temporary
//' # directory. The data are deleted from the disk as soon as the object
//' # is destroyed.
//' sim <- TissueSimulation("test")
//'
//' # the name of the simulation is "test"
//' sim$get_name()
//'
//' # however no directory "test" has been created in the working directory
//' "test" %in% list.files()
//'
//' # By using the optional parameter `save_directory`, we force the
//' # simulation to save its data in a local directory whose name
//' # is the name of the simulation, i.e., "test". This data will be
//' # preserved when the simulation object will be destroyed.
//' sim <- TissueSimulation("test", save_directory = TRUE)
//'
//' # the directory "test" exists and contains a binary dump of
//' # the simulation
//' "test" %in% list.files()
//'
//' # the directory persists even after the object destruction
//' rm(sim)
//' "test" %in% list.files()
//'
//' # let us manually delete the "test" directory
//' unlink("test", recursive = TRUE)
//'
//' # the name parameter is optional
//' sim <- TissueSimulation(save_directory = TRUE)
//'
//' # the name of the simulation is `ProCESS_<YY><MM><DD>_<HH><MM><SS>`
//' sim$get_name()
//'
//' # the simulation dump have been saved in a directory named
//' # after the simulation name
//' list.files(pattern = "^ProCESS_")
//'
//' # let us remove the object and manually delete the simulation
//' # directory
//' rm(sim)
//' unlink(list.files(pattern = "^ProCESS_"), recursive = TRUE)
//'
//' # users can provide a random seed to the simulation...
//' sim <- TissueSimulation(seed = 13)
//'
//' # ..., specify the size of the simulated space by using the
//' # optional parameters `width` and `height`, or...
//' sim <- TissueSimulation(width = 1200, height = 900)
//' sim$get_tissue_size()
//'
//' # ... build a simulation, add its species, and set their rates
//' # by passing a data frame specifying the event rates
//' df_rates <- data.frame(
//'   mutant = c("A", "A", "B"),
//'   event = c("duplication", "death", "duplication"),
//'   rate = c(0.3, 0.2, 0.2)
//' )
//' df_rates
//'
//' sim <- TissueSimulation(rates = df_rates)
//' sim
//'
//' # if epigenetic states are needed, the data frame must also contain
//' # the columns `epistate` and `first.child.epistate`
//' df_rates <- data.frame(
//'   mutant = c("A", "A", "A", "A", "A", "B", "B"),
//'   epistate = c("E1", "E1", "E1", "E1", "E2", "E1", "E1"),
//'   event = c("duplication", "death", "switch", "switch",
//'             "duplication", "duplication", "switch"),
//'   first.child.epistate = c("E1", NA, "E2", "E3", "E2", "E1", "E3"),
//'   rate = c(0.3, 0.2, 0.01, 0.04, 0.2, 0.2, 0.01)
//' )
//' df_rates
//'
//' sim <- TissueSimulation(rates = df_rates)
//' sim
//' @seealso <code>[TissueSimulation](TissueSimulation_class.md)</code>
    function("TissueSimulation", &TissueSimulation::build_simulation,
             List::create(_["name"] = R_NilValue, _["width"] = 1000, _["height"] = 1000,
                          _["save_directory"] = false, _["rates"] = R_NilValue,
                          _["epigenetic_states"] = R_NilValue,
                          _["seed"] = R_NilValue),
             "Create a tissue simulation");

//' @name SampleForest
//' @title The sample cell ancestor forest
//' @description The ancestor forest of sampled cells
//' @details This class represents the forest of the ancestors of the
//'   cells sampled during the computation. The leaves of
//'   this forest are the sampled cells.
//'
//'   The objects of this class provide the following methods and properties:
//'   - <code>[get_coalescent_cells()](SampleForest-cash-get_coalescent_cells.md)</code>
//'     returns the most recent common ancestors of the sampled cells.
//'   - <code>[get_node()](SampleForest-cash-get_node.md)</code>
//'     returns an object of type <code>[SampleForestNode]</code>.
//'   - <code>[get_nodes()](SampleForest-cash-get_nodes.md)</code>
//'     returns information about the nodes in the forest.
//'   - <code>[get_samples_info()](SampleForest-cash-get_samples_info.md)</code>
//'     returns information about the samples generating the forest.
//'   - <code>[get_species_info()](SampleForest-cash-get_species_info.md)</code>
//'     returns information about the simulated species.
//'   - <code>[get_sticks()](SampleForest-cash-get_sticks.md)</code>
//'     computes the forest sticks.
//'   - <code>[get_subforest_for()](SampleForest-cash-get_subforest_for.md)</code>
//'     extracts a sub-forest for some of the samples.
//'   - <code>[represents_cell()](SampleForest-cash-represents_cell.md)</code>
//'     tests whether a cell having a given identifier is represented by the
//'     forest.
//'   - <code>[save()](SampleForest-cash-save.md)</code>
//'     saves the forest.
//' @seealso <code>[PhylogeneticForest]</code>
    class_<SampleForest>("SampleForest")
        REGISTER_FOREST_COMMON_FIELD(SampleForest);

//' @name SampleForestNode
//' @title The node of a sample forest
//' @description A class representing the nodes of a sample forest.
//' @details This class represents the nodes of a sample forest. It does not
//'   have a user constructor. Its objects are produced by [get_node_tour()] and
//'   <code>[SampleForest$get_node()]</code>.
//'
//'   The objects of this class provide the following methods and properties:
//'   - <code>[cell_id](SampleForestNode-cash-cell_id.md)</code>
//'     represents the identifier of the associated cell.
//'   - <code>[parent](SampleForestNode-cash-parent.md)</code>
//'     represents the node's parent.
//'   - <code>[children](SampleForestNode-cash-children.md)</code>
//'     represents a list of the node's children.
//'   - <code>[is_root](SampleForestNode-cash-is_root.md)</code>
//'     is a Boolean flag that is `TRUE` if and only if the node is a forest
//'     root.
//'   - <code>[is_leaf](SampleForestNode-cash-is_leaf.md)</code>
//'     is a Boolean flag that is `TRUE` if and only if the node is a forest
//'     leaf.
//'   - <code>[sample_name](SampleForestNode-cash-sample_name.md)</code>
//'     is the name of the sample that collected the associated cell.
//'   - <code>[birth_time](SampleForestNode-cash-birth_time.md)</code>
//'     is the birth time of the cell associated to the node.
//'   - <code>[death_time](SampleForestNode-cash-death_time.md)</code>
//'     is the death time of the cell associated to the node.
//'   - <code>[life_span](SampleForestNode-cash-life_span.md)</code>
//'     is the life span of the cell associated to the node.
//'   - <code>[species_name](SampleForestNode-cash-species_name.md)</code>
//'     is the name of the associated cell's species.
//'   - <code>[epistate_name](SampleForestNode-cash-epistate_name.md)</code>
//'     is the name of the associated cell's epigenetic state.
//'   - <code>[mutant_name](SampleForestNode-cash-mutant_name.md)</code>
//'     is the name of the associated cell's mutant.
//'
//' @seealso [get_node_tour()], <code>[SampleForestNodeTour]</code>,
//'   <code>[PhylogeneticForestNode]</code>,
//'   `vignette("node_tour")`
   class_<SampleForestNode>("SampleForestNode")
        REGISTER_NODE_COMMON_FIELDS(SampleForestNode);

//' @name SampleForestNodeTour
//' @title An iterator class over sample forest nodes
//' @description Iterators over sample forest nodes.
//' @details This class represents iterators over sample forest nodes.
//'   The objects of this class are built by [get_node_tour()] and provide
//'   the following methods and properties:
//'   - <code>[node](SampleForestNodeTour-cash-node.md)</code>
//'     is an object of the class <code>[SampleForestNode]</code> and
//'     represents the node pointed by the iterator.
//'   - <code>[label](SampleForestNodeTour-cash-label.md)</code> (OPTIONAL)
//'     represents the label of the of the node pointed by the
//'     iterator. The presence of this field depends on the type of the
//'     [get_node_tour()]'s parameters used to create the tour object.
//'   - <code>[step()](SampleForestNodeTour-cash-step.md)</code>
//'     moves the iterator to the next node in the tour.
//'   - <code>[done](SampleForestNodeTour-cash-done.md)</code>
//'     is a Boolean flag that is set to `TRUE` only when the tour ended.
//'
//' @seealso [get_node_tour()], <code>[SampleForestNode]</code>,
//'   <code>[PhylogeneticForestNodeTour]</code>,
//'   `vignette("node_tour")`
    class_<SampleForestNodeTour>("SampleForestNodeTour")
        REGISTER_TOUR_COMMON_FIELDS(SampleForestNodeTour);

//' @name get_node_tour
//' @title Labelling forest nodes
//' @description This method generates a <code>[SampleForestNodeTour]</code>
//' @usage get_node_tour(forest, labelling_functor, init_value, only_leaves,
//'                       with_genomes)
//' @param forest Either a <code>[SampleForest]</code> or a
//'  <code>[PhylogeneticForest]</code> object.
//' @param labelling_functor Depending on the type of `forest`, a function
//'  of the type `label_type (*)(label_type, SampleForestNode)`, when
//'  `forest` is a <code>[SampleForest]</code> object, or
//'  `label_type (*)(label_type, PhylogeneticForestNode)`, when `forest`
//'  is is a <code>[PhylogeneticForest]</code> object (default: `NULL`).
//' @param init_value The initial value of the labelling process
//'  (default: `NULL`).
//' @param only_leaves A Boolean flag to iterate exclusively over the
//'   forest leaves (default: `FALSE`).
//' @param with_genomes A Boolean flag to also generate the corresponding
//'   cell genomes. It can be set to TRUE exclusively when `forest` is
//'   a <code>[PhylogeneticForest]</code> object. (default: `FALSE`)
//' @return Either a <code>[SampleForestNodeTour]</code> or a
//'   <code>[PhylogeneticForestNodeTour]</code> iterating over `forest`'s
//'   nodes. When `only_leaves` is set to `TRUE`, the returned tour
//'   iterates over `forest`'s leaves.
//' @examples
//' # use a sample forest example
//' forest <- example("SampleForest")
//'
//' # get a tour of the forest nodes
//' node_tour <- get_node_tour(forest)
//'
//' # the first node in the tour
//' node_tour$node
//'
//' # the first node is not a leaf
//' node_tour$node$is_leaf
//'
//' # move to the next node
//' node_tour$step()
//'
//' # not even the second node is a leaf
//' node_tour$node$is_leaf
//'
//' @examples
//' # use a sample forest example
//' forest <- example("SampleForest")
//'
//' # get a tour of the forest leaves
//' node_tour <- get_node_tour(forest, only_leaves = TRUE)
//'
//' # the first node in the tour is a leaf
//' node_tour$node
//' node_tour$node$is_leaf
//'
//' # move to the next node
//' node_tour$step()
//'
//' # the second node too
//' node_tour$node$is_leaf
//'
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get a tour of the forest nodes with their genomes
//' node_tour <- get_node_tour(forest, with_genomes = TRUE)
//'
//' # the first node in the tour is a leaf
//' node_tour$node
//'
//' # the forest was build using the setup "demo". Thus, the
//' # node genomes have only chromosome 22
//' node_tour$genome
//'
//' @examples
//' # use a sample forest example
//' forest <- example("SampleForest")
//'
//' # we define a function to collect the tour labels
//' collect_labels <- function(tour) {
//'   total <- NULL
//'
//'   # `tour$done` is `TRUE` iff the tour ended
//'   while (!tour$done) {
//'     if (is.null(total)) {
//'       # `tour$label` is the node label
//'       total <- list(tour$label)
//'     } else {
//'       total <- append(total, tour$label)
//'     }
//'
//'     # `tour$step()` advances to the next node
//'     # in the tour
//'     tour$step()
//'   }
//'
//'   total
//' }
//'
//' print("Functor 1")
//' # a labelling functor
//' labelling_functor1 <- function(label, node) {
//'   # the nodes are labelled by the identifiers of the associated cells
//'   node$cell_id
//' }
//'
//' # since `labelling_functor1` does not use `label`, we can omit the
//' # parameter `init_value`
//' tour <- get_node_tour(forest, labelling_functor1)
//'
//' print("Functor 1 - All nodes")
//' print(collect_labels(tour)[1:5])
//'
//' # since `labelling_functor1` does not use `label`, we can omit the
//' # parameter `init_value`
//' tour <- get_node_tour(forest, labelling_functor1, only_leaves = TRUE)
//'
//' print("Functor 1 - Only leaves")
//' print(collect_labels(tour)[1:5])
//'
//' labelling_functor2 <- function(label, node) {
//'   # the nodes are labelled by their visiting order
//'   label + 1
//' }
//'
//' # `labelling_functor2` uses `label` and we must specify the
//' # parameter `init_value`
//' tour <- get_node_tour(forest, labelling_functor2,
//'                        init_value = 0, only_leaves = TRUE)
//'
//' print("Functor 2 - Only leaves")
//' print(collect_labels(tour)[1:5])
//'
//' a <- 3.14
//'
//' # this functor uses a global variable to compute the
//' # labels
//' labelling_functor3 <- function(label, node) {
//'   # the nodes are labelled by their visiting order multiplied
//'   # by the value in the global variable `a`
//'   label + a
//' }
//'
//' tour <- get_node_tour(forest, labelling_functor3,
//'                        init_value = 0, only_leaves = TRUE)
//'
//' print("Functor 3 - Only leaves")
//' print(collect_labels(tour)[1:5])
//'
//' set.seed(0)
//'
//' # this functor uses a random function
//' labelling_functor4 <- function(label, node) {
//'   # the nodes are randomly labelled
//'   label + sample(-100:100, 1)
//' }
//'
//' tour <- get_node_tour(forest, labelling_functor4,
//'                        init_value = 0, only_leaves = TRUE)
//'
//' print("Functor 4 - Only leaves")
//' print(collect_labels(tour)[1:5])
//' @seealso <code>[SampleForestNode]</code>,
//'    <code>[SampleForestNodeTour]</code>,
//'    <code>[PhylogeneticForestNode]</code>,
//'    <code>[PhylogeneticForestNodeTour]</code>,
//'    `vignette("node_tour")`
    function("get_node_tour", &get_node_tour,
            List::create(
                 _["forest"], _["labelling_functor"] = R_NilValue,
                 _["init_value"] = R_NilValue, _["only_leaves"] = false,
                 _["with_genomes"] = false),
            "Get a labelling tour");

//' @name load_sample_forest
//' @title Loading sample forests
//' @description This method loads a sample forest from a file.
//' @usage load_sample_forest(filename, quiet)
//' @param filename The path of the file from which the samples
//'   forest must be loaded.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @return The loaded sample forest
//' @seealso [SampleForest$save()], [load_phylogenetic_forest()]
//'   [load_forest()]
    function("load_sample_forest",
             (SampleForest (*)(const std::string &,
                                     const bool))&SampleForest::load,
             List::create(_["filename"] = "", _["quiet"] = false),
             "Load a sample forest");

//' @name load_forest
//' @title Loading forests
//' @description This method loads a forest from a file.
//' @details This method loads a forest, being either a sample forest or a
//'   phylogenetic forest, from a file.
//' @usage load_forest(filename, quiet)
//' @param filename The path of the file from which the forest must be
//'   loaded.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: `FALSE`).
//' @return The loaded forest
//' @seealso [SampleForest$save()], [PhylogeneticForest$save()],
//'   [load_sample_forest()], [load_phylogenetic_forest()]
    function("load_forest",
             (SEXP (*)(const std::string &,
                       const bool))&ForestCore::load_forest,
             List::create(_["filename"] = "", _["quiet"] = false),
             "Load a forest");
}
