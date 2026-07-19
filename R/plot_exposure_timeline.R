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

get_exposure_ends <- function(phylo_forest) {
  exposure <- phylo_forest$get_exposures()

  time_points <- exposure%>% dplyr::pull(time) %>% unique
  signatures <- exposure %>% dplyr::pull(signature) %>% unique
  for (t in time_points) {
    for (signature in signatures) {
      if (nrow(exposure %>%
                 dplyr::filter(.data$time == t,
                               .data$signature == signature)) == 0) {
        exposure[nrow(exposure)+1,] <- c(as.numeric(t), signature,
                                         as.numeric(0), gsub("[0-9]*$","",
                                                             signature))
      }
    }
  }

  last_time <- max(phylo_forest$get_samples_info()["time"])
  time_points <- sort(unique(append(time_points, last_time)))
  end_time <- apply(exposure, 1, function(x) {
    time_points[which(round(time_points, 2) == round(as.numeric(x[1]), 2)) + 1]
  })
  end_time[is.na(end_time)] <- as.numeric(last_time)

  exposure$end_time <- end_time

  exposure <- exposure[order(exposure$time, exposure$signature), ]

  exposure[, c(1, 3)] <- sapply(exposure[, c(1, 3)], as.numeric)

  return(exposure)
}

#' Plot the signature exposure timeline of a phylogenetic forest
#'
#' @description
#' Plots the signatures exposure changes along a phylogenetic
#' forest.
#'
#' @param phylogenetic_forest A phylogenetic forest.
#' @param linewidth The width of the lines in the plot.
#' @param emphasize_switches A Boolean flag to emphasize the
#'   exposure switches.
#' @param pal_name The name of a `RColorBrewer` palette.
#'
#' @return A `ggplot2` plot.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # plotting the phylogenetic forest
#' plot_exposure_timeline(forest)
#'
#' # plotting the phylogenetic forest emphasizing the exposure switches
#' plot_exposure_timeline(forest, emphasize_switches = TRUE)
#' @export
plot_exposure_timeline <- function(phylogenetic_forest, linewidth = 0.8,
                                   emphasize_switches = FALSE,
                                   pal_name = "Set3") {
  stopifnot(inherits(phylogenetic_forest, "Rcpp_PhylogeneticForest"))

  exposure_df <- get_exposure_ends(phylogenetic_forest)

  signames <- exposure_df %>% dplyr::pull(signature) %>% unique
  colors <- get_colors_for(signames, pal_name = pal_name)

  line_shift <- 0.005 * linewidth

  exposure_df <- exposure_df %>% dplyr::group_by(time,exposure) %>%
    dplyr::mutate(line_pos = .data$exposure +
                    line_shift * (dplyr::row_number() - 1))

  p <- ggplot2::ggplot(data = exposure_df,
                       ggplot2::aes(x = time, y = .data$line_pos,
                                    group = signature, color = signature)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$end_time,
                                       yend = .data$line_pos),
                          linewidth = linewidth) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Time Point", y = "Signature Exposure") +
    ggplot2::scale_colour_manual(name = "Signatures", values = colors) +
    ggplot2::theme_minimal()  # Apply a minimal theme

  if (emphasize_switches) {
    switching_times <- exposure_df["time"] %>% dplyr::filter(time != 0)
    for (switching_time in switching_times) {
      p <- p +
        ggplot2::geom_vline(xintercept = switching_time, linetype = "dotted")
    }
  }

  return(p)
}
