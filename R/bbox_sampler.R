## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2026 - Giulio Caravagna <gcaravagna@units.it>
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

#' Bounding box sampler
#'
#' @description
#' This function searches for a tissue rectangle that contains a given
#' number of cells of a certain type. The search is performed by sampling
#' the tissue uniformly until the constraint on the number of cells
#' contained in the rectangle is satisfied or the sampling has been
#' repeated for a specified number of times. When the constraint cannot be
#' satisfied, the rectangle having the maximal number of cells of the
#' specified type among those sampled is returned.
#'
#' @param simulation A simulation object.
#' @param which The species name.
#' @param n The desired number of cells from species `which`.
#' @param n_w Width of the box.
#' @param n_h Height of the box.
#' @param nattempts The maximum number of attempts.
#'
#' @return coordinates for a bounding box.
#' @export
#'
#' @examples
#' sim <- TissueSimulation()
#' sim$add_mutant("A", c(duplication = 0.08, death = 0.01))
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_size("A", 25000)
#' bbox <- bbox_sampler(sim, "A", n = 2500, n_w = 50, n_h = 50)
#' sim$sample_cells("A", bbox$p, bbox$q)
#' plot_tissue(sim)
bbox_sampler <- function(simulation, which, n, n_w, n_h, nattempts = 100) {

  sp <- cli::make_spinner()

  # where are cells of species which
  locations <- simulation$get_cells() %>%
    add_species_col() %>%
    dplyr::filter(.data$species == which) %>%
    dplyr::select(.data$position_x, .data$position_y)

  bof_p <- bof_q <- NULL
  bof_counts <- -1

  repeat {
    sp$spin()

    location_id <- sample(seq_len(nrow(locations)), 1)

    p <- locations[location_id, ] %>% as.numeric()
    q <- p + c(n_w, n_h)

    nc <- simulation$get_cells(p, q) %>%
      add_species_col() %>%
      dplyr::filter(.data$species == which) %>%
      nrow()

    if (nc > bof_counts) {
      bof_p <- p
      bof_q <- q
      bof_counts <- nc
    }

    if (nc >= n || nattempts == 0) break
    nattempts <- nattempts - 1
  }

  sp$finish()

  return(list(p = bof_p, q = bof_q))
}
