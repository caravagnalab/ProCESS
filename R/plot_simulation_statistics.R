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

#' Plot the number of stochastic events in the simulation
#'
#' @description
#' A pie chart with events split by type, mutant and epigentic state
#' where they occurred. It also provides annotations for the simulation
#' information.
#'
#' @param simulation A simulation.
#'
#' @return A ggplot plot.
#' @export
#' @examples
#' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
#' sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
#'                                    E2 = 0.01),
#'                          E2 = list(duplication = 0.08, death = 0.01,
#'                                    E1 = 0.02)))
#' sim$place_cell("A[E1]", 500, 500)
#' sim$run_up_to_time(60)
#' plot_firings(sim)
plot_firings <- function(simulation) {
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))

  firings <- simulation$get_firings() %>%
    add_species_col()

  plot <- ggplot2::ggplot(firings) +
    ggplot2::geom_bar(stat = "identity",
                      ggplot2::aes(x = "", y = .data$fired,
                                   fill = .data$event))
  if ("epistate" %in% colnames(firings)) {
    plot <- plot + ggplot2::facet_grid(.data$mutant ~ .data$epistate)
  } else {
    plot <- plot + ggplot2::facet_grid(.data$mutant)
  }

  plot + ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Event",
      x = NULL,
      y = NULL
    ) +
    my_theme() +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme(legend.position = "bottom")
}