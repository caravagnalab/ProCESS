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
                              max_mutants = 3) {

  if (is.null(.pkg_env$mutant_color_map)) {
    max_mutants <- max(max_mutants, length(mutants))
    .pkg_env$mutant_color_map <- RColorBrewer::brewer.pal(max_mutants, pal_name)
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

get_species_colors <- function(data, pal_name = "Dark2",
                               max_mutants = 4) {
  if (inherits(data, "Rcpp_TissueSimulation")) {
    data <- data$get_counts() %>%
      dplyr::select(-.data$counts, -.data$overall)
  } else if (inherits(data, "Rcpp_SampleForest")
             ||inherits(data, "Rcpp_PhylogeneticForest")) {
    data <- data$get_species_info()
  }

  mutants <- data %>%
    dplyr::pull(.data$mutant) %>% unique()

  mutant_color_map <- get_mutant_colors(mutants, pal_name = pal_name,
                                        max_mutants = max_mutants)

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