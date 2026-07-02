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

my_theme <- function() {
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

get_colors_for <- function(values, pal_name = "Dark2") {
  # Unpaired
  colors <- NULL
  if (length(values) > 0) {
    num_of_values <- values %>% length()

    if (num_of_values < 3) {
      colors <- RColorBrewer::brewer.pal(3, pal_name)

      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }

    names(colors) <- values
  }

  return(colors)
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

  new_data %>%
    dplyr::arrange(.data$mutant)
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

get_species_colors <- function(data) {
  if ("epistate" %in% colnames(data)) {
    name_species <- data %>%
      dplyr::mutate(
        species = species_name(.data$mutant, .data$epistate)
      ) %>%
      dplyr::arrange(.data$mutant)

    get_colors_for(name_species[, "species"], "Paired")
  } else {
    name_species <- data %>%
      dplyr::arrange(.data$mutant)

    get_colors_for(name_species[, "mutant"], "Dark2")
  }
}

validate_chromosomes <- function(seq_res, chromosomes) {
  seq_res_chrs <- (seq_res["chr"] %>% unique())[,1]
  if (is.null(chromosomes)) {
    chromosomes <- seq_res_chrs
  } else {
    unknown_chrs <- dplyr::setdiff(chromosomes, seq_res_chrs)

    if (length(unknown_chrs)>0) {
        unknown_chrs_str <- paste0(unknown_chrs, collapse = ", ")
        if (length(unknown_chrs)>1) {
            msg <- paste0("The chromosomes ", unknown_chrs_str, " are")
        } else {
            msg <- paste0("The chromosome ", unknown_chrs_str, " is")
        }
        stop(paste0(msg, " not present in the sequence reference data."))
    }
  }

  return(chromosomes)
}