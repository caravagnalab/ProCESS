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

my_theme <- function() {
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

get_colors_for <- function(values, pal_name = "Dark2",
                           min_values = 4) {
  # Unpaired
  colors <- NULL
  if (length(values) > 0) {
    num_of_values <- values %>% length()

    if (num_of_values < min_values) {
      colors <- RColorBrewer::brewer.pal(min_values, pal_name)

      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }

    names(colors) <- values
  }

  colors
}

get_mutant_colors <- function(mutants, pal_name = "Dark2",
                              min_mutants = 3) {

  if (is.null(.pkg_env$mutant_color_map)) {
    min_mutants <- max(min_mutants, length(mutants))
    .pkg_env$mutant_color_map <- RColorBrewer::brewer.pal(min_mutants, pal_name)
  } else {
    if (length(mutants) > length(.pkg_env$mutant_color_map)) {
      warning(paste0("The number of mutants changed since",
                     " the last plot. The color map has changed."))

      .pkg_env$mutant_color_map <- RColorBrewer::brewer.pal(length(mutants),
                                                            pal_name)
    }
  }

  colors <- .pkg_env$mutant_color_map[seq_along(mutants)]

  names(colors) <- mutants

  colors
}

species_name <- function(mutant, epistate) {
  dplyr::if_else(is.na(epistate),
                 epistate,
                 paste0(mutant, "[", epistate, "]"))
}

add_species_col <- function(data, col_name = "species") {
  if ("epistate" %in% colnames(data)) {
    new_data <- data %>% dplyr::mutate(
      !!dplyr::sym(col_name) := species_name(.data$mutant, .data$epistate)
    )
  } else {
    new_data <- data %>% dplyr::mutate(
      !!dplyr::sym(col_name) := .data$mutant
    )
  }

  new_data
}

get_species <- function(simulation) {
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))

  with_epigenetics <- nrow(simulation$get_epigenetic_states()) > 0

  if (with_epigenetics) {
    dplyr::cross_join(simulation$get_mutants(),
                      simulation$get_epigenetic_states())
  } else {
    simulation$get_mutants()
  }
}

get_epistate_shade_map <- function(epistates, min_light_value = 0.3) {
  if (!is.null(.pkg_env$epistate_shade_map)) {
    if (length(.pkg_env$epistate_shade_map) < length(epistates)) {
      warning(paste0("The number of epigenetic states changed since",
                     " the last plot. The color map has changed."))

      .pkg_env$epistate_shade_map <- seq(min_light_value, 0,
                                         length.out = length(epistates))

      names(.pkg_env$epistate_shade_map) <- epistates
    }
  } else {
    .pkg_env$epistate_shade_map <- seq(min_light_value, 0,
                                       length.out = length(epistates))

    names(.pkg_env$epistate_shade_map) <- epistates
  }

  .pkg_env$epistate_shade_map
}

#' Get species color maps
#'
#' @description
#' This function returns a species color map
#' @details
#' This function returns a color maps for the species represented
#' in a tissue, in a forest, or in a data frame containing
#' mutants and, possible, epigenetic states.
#' @param object The tissue simulation, the forest, or the data
#'   frame whose species color map is required.
#' @param pal_name The `RColorBrewer` palette name used to generate
#'   the species color map (default: `Dark2`).
#' @param min_mutants The minimum number of colors generated.
#' @return A named list whose name are the names of the species
#'   in `object` and whose values are the associated colors.
#' @examples
#' # get an example of `SampleForest` object
#' forest <- example("SampleForest")
#'
#' # get the species info
#' forest$get_species_info()
#'
#' # get the color map for the forest species
#' get_species_colors(forest)
#' @export
#'
get_species_colors <- function(object, pal_name = "Dark2",
                               min_mutants = 4) {
  if (inherits(object, "Rcpp_TissueSimulation")) {
    data <- object$get_counts() %>%
      dplyr::select(-.data$counts, -.data$overall)
  } else if (inherits(object, "Rcpp_SampleForest")
             ||inherits(object, "Rcpp_PhylogeneticForest")) {
    data <- object$get_species_info()
  } else {
    if (!is.data.frame(object)) {
      stop(paste("The 1st parameter must be a `TissueSimulation`, a",
                 "forest, or a data frame."))
    }
    if (!"mutant" %in% names(object)) {
      stop("The 1st parameter does not contain the column \"mutant\".")
    }
    data <- object
  }

  mutants <- data %>%
    dplyr::pull(.data$mutant) %>% unique()

  mutant_color_map <- get_mutant_colors(mutants, pal_name = pal_name,
                                        min_mutants = min_mutants)

  if ("epistate" %in% colnames(data)) {

    epistates <- data %>%
      dplyr::pull(.data$epistate) %>% unique()

    epistate_shade_map <- get_epistate_shade_map(epistates)

    species_names <- c()
    species_color_map <- c()

    for (i in seq_along(mutants)) {
      for (j in seq_along(epistate_shade_map)) {
        color <- colorspace::lighten(mutant_color_map[[mutants[i]]],
                                     amount = epistate_shade_map[j])
        species_color_map <- c(species_color_map, color)
        species_names <- c(species_names,
                           species_name(mutants[i], epistates[j]))
      }
    }

    names(species_color_map) <- species_names

    species_color_map
  } else {
    mutant_color_map
  }
}

validate_chromosomes <- function(seq_res, chromosomes) {
  seq_res_chrs <- (seq_res["chr"] %>% unique())[, 1]
  if (is.null(chromosomes)) {
    chromosomes <- seq_res_chrs
  } else {
    unknown_chrs <- dplyr::setdiff(chromosomes, seq_res_chrs)

    if (length(unknown_chrs) > 0) {
      unknown_chrs_str <- paste0(unknown_chrs, collapse = ", ")
      if (length(unknown_chrs) > 1) {
        msg <- paste0("The chromosomes ", unknown_chrs_str, " are")
      } else {
        msg <- paste0("The chromosome ", unknown_chrs_str, " is")
      }
      stop(paste0(msg, " not present in the sequence reference data."))
    }
  }

  chromosomes
}