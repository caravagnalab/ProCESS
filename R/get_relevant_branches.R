## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2025 - Giulio Caravagna <gcaravagna@units.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Get the data frame of the relevant branches
#'
#' @description
#' Get subset of forest nodes corresponding to phylogenetic branches containing
#'  relevant biological events (drivers)
#'
#' @param forest A forest
#'
#' @return A data frame reporting all the forest nodes corresponding to
#'  phylogenetic branches containing relevant biological events (drivers).
#' @export
#'
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- TissueSimulation()
#'
#' sim$add_mutant("A", c(duplication = 1))
#' sim$place_cell("A", 500, 500)
#'
#' sim$run_up_to_size("A", 1e4)
#'
#' sim$add_mutant("B", c(duplication = 3.5))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B", 1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_sample_forest()
#' get_relevant_branches(forest)

get_relevant_branches <- function(forest) {
  if (!inherits(forest, "Rcpp_SampleForest") &&
      !inherits(forest, "Rcpp_PhylogeneticForest")) {
    stop("The parameter must be a forest.")
  }

  sticks <- forest$get_sticks()

  if (length(sticks) > 0) {
    nodes <- forest$get_nodes() %>%
      dplyr::select(-.data$depth, -.data$birth_time, -.data$sample)

    events_table <- lapply(seq_along(sticks), function(i) {

      first_in_stick <- nodes %>%
        dplyr::filter(.data$cell_id == sticks[[i]][1]) %>%
        add_species_col()

      label <- first_in_stick %>% dplyr::pull(.data$species)

      data.frame(cell_id = sort(sticks[[i]][2:length(sticks[[i]])]),
                 label = label)
    }) %>% dplyr::bind_rows() %>% unique()

    forest$get_nodes() %>%
      dplyr::full_join(events_table, by = "cell_id") %>%
      dplyr::mutate(label = ifelse(is.na(.data$label),
                                   "Subclonal", .data$label)) %>%
      dplyr::mutate(label = ifelse(.data$cell_id == 1,
                                   "Truncal", .data$label))

  } else {
    NULL
  }

}
