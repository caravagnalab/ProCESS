## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2026 - Giulio Caravagna <gcaravagna@units.it>
##                           Alberto Casagrande <alberto.casagrande@uniud.it>
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

#' Plot the number of cells along a simulation
#'
#' @description
#' A pie chart with population counts, split by mutants and epistate. It
#' also provides annotations for the simulation data.
#'
#' @param simulation A simulation.
#' @param color_map A named vector representing the simulation species color
#'   map (optional).
#' @return A ggplot plot.
#' @export
#'
#' @examples
#' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
#' sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
#'                                    E2 = 0.01),
#'                          E2 = list(duplication = 0.08, death = 0.01,
#'                                    E1 = 0.02)))
#' sim$history_delta <- 1
#' sim$place_cell("A[E1]", 500, 500)
#' sim$run_up_to_time(60)
#'
#' plot_timeseries(sim)
#'
#' # define a custom color map
#' color_map <- c("#B2DF8A", "#E31A1C")
#' names(color_map) <- c("A[E1]", "A[E2]")
#'
#' plot_timeseries(sim, color_map = color_map)
plot_timeseries <- function(simulation, color_map = NULL) {
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))

  counts <- simulation$get_count_history() %>%
    add_species_col()
  time <- simulation$get_clock() %>% round(digits = 3)

  if (is.null(color_map)) {
    color_map <- get_species_colors(simulation)
  }

  counts$species <- factor(counts$species,
                           levels = unique(names(color_map)))

  ggplot2::ggplot(counts) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = .data$count,
                                    color = .data$species)) +
    ggplot2::geom_point(ggplot2::aes(x = time, y = .data$count,
                                     color = .data$species)) +
    ggplot2::labs(
      x = "Time", y = "Cell count",
      color = "Species"
    ) +
    my_theme() +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::theme(legend.position = "bottom")
}