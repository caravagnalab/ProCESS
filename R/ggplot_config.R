## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023 - Giulio Caravagna <gcaravagna@units.it>
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

get_species_colors <- function(species) {

  paired_species <- species %>%
    dplyr::filter(!is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$mutant, .data$epistate)
    ) %>%
    dplyr::arrange(.data$mutant)

  unpaired_species <- species %>%
    dplyr::filter(is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$mutant, .data$epistate)
    ) %>%
    dplyr::arrange(.data$mutant)

  # Unpaired
  unp_colour <- NULL
  if (nrow(unpaired_species) > 0) {
    unp_colour <- get_colors_for(unpaired_species[, "mutant"], "Dark2")
  }

  # Paired
  pai_colour <- NULL
  if (nrow(paired_species) > 0) {
    pai_colour <- get_colors_for(paired_species[, "species"], "Paired")
  }

  return(c(pai_colour, unp_colour))
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