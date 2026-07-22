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

#' Plot a tissue
#'
#' @description
#' Plots cells distribution over a tissue highlighting species by color.
#' To facilitate the plot and avoid excessive number of cells, for instance,
#' when a simulation deals with millions of cells, the plot draws a
#' hexagonal heatmap of 2D bins.
#'
#' @param simulation A simulation object.
#' @param num_of_bins The number of bins (default: 100).
#' @param before_sample A sample name. When provided, this function represents
#'   the tissue as appeared just before the specified sampling. The parameters
#'   `before_sample` and `at_sample` are mutually exclusive (optional).
#' @param at_sample A sample name. When provided, this function represents
#'   the tissue as appeared when the specified sampling was about to be
#'   collected. The parameters `before_sample` and `at_sample` are mutually
#'   exclusive (optional).
#' @param plot_next_sample_regions A Boolean value. When `before_sample` is
#'   set and `plot_next_sample_regions` is set to be `TRUE`, this function
#'   plots the regions of the samples collected at the same simulated time
#'   of the specified sample. When, instead, `at_sample` is set and
#'   `plot_next_sample_regions` is set to be `TRUE`, the function plots the
#'   regions of the samples collected at the same simulated time of the
#'   specified sample, but not before the specified sample (default: `FALSE`).
#' @param plot_sample_region A Boolean value. When either `at_sample` or
#'   `before_sample` are set and `plot_sample_region` is set to be `TRUE`,
#'   the function also plots the region of the specified sample
#'   (default: `TRUE`).
#' @param color_map A named vector representing the simulation species color
#'   map (optional).
#' @param list_all_species A Boolean flag to show all species in
#'   the legend (default: `FALSE`).
#' @param highlight_function A function that takes as the input each row of
#'   the data frame returned by <code>[TissueSimulation$get_cells()]</code>
#'   and returns a Boolean value. When the function returns `FALSE` the
#'   corresponding cells is plotted in grey. When the parameter is set to
#'   `NULL`, all tumour simulation cells are colored (default: `NULL`).
#' @return An editable ggplot plot.
#' @examples
#' # set the seed
#' set.seed(0)
#'
#' # build a tissue simulation
#' sim <- TissueSimulation(width = 600, height = 600)
#'
#' # add the mutant A
#' sim$add_mutant("A", c(duplication = 0.12, death = 0.05))
#'
#' # place a cell in the tissue and simulate it until 10 cells
#' sim$place_cell("A", 300, 300)
#' sim$run_up_to_size("A", 10)
#'
#' # add the mutant B and let mutate a border cell of A in B
#' sim$add_mutant("B", c(duplication = 0.145, death = 0.06))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
#'
#' # simulate the tissue up to 30 cells in B
#' sim$run_up_to_size("B", 30)
#'
#' # add the third mutant and let one cell of A mutate into C
#' sim$add_mutant("C", c(duplication = 0.15, death = 0.06))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "C")
#'
#' # simulate the tissue until C consists of 25000 cells
#' sim$run_up_to_size("C", 25000)
#'
#' # collect the sample "S1"
#' sim$sample_cells("S1", c(145, 230), c(215, 300))
#'
#' # let the simulation reach 25000 cells in C again
#' sim$run_up_to_size("C", 25000)
#'
#' # collect two samples
#' sim$sample_cells("S2", c(350, 300), c(420, 370))
#' sim$sample_cells("S3", c(200, 350), c(270, 420))
#'
#' # add a further mutant and derive it from B
#' sim$add_mutant(name = "D", c(duplication = 0.8, death = 0.01))
#' sim$mutate_progeny(sim$choose_border_cell_in("B"), "D")
#'
#' # let the tumour evolve until the mutant C and D cumulatively
#' # consist of 10000 cells
#' sim$run_until(sim$var("C") + sim$var("D") == 1e5)
#'
#' # plot the tissue in the current status
#' plot_tissue(sim)
#'
#' # plot the tissue as it was when "S3" was about to be sampled
#' plot_tissue(sim, at_sample = "S3")
#'
#' # plot the tissue as it was when "S3" was about to be sampled and
#' # highlight the regions of the samples collected at the same simulated
#' # time, but not before
#' # it, i.e., "S3"
#' plot_tissue(sim, at_sample="S3",
#'             plot_next_sample_regions = TRUE)
#'
#' # plot the tissue as it was just before sampling "S3"
#' plot_tissue(sim, before_sample="S3")
#'
#' # plot the tissue as it was just before sampling "S3" and highlight
#' # the regions of the samples collected at the same simulated time,
#' # i.e., "S2" and "S3"
#' plot_tissue(sim, before_sample="S3",
#'             plot_next_sample_regions = TRUE)
#'
#' # define a custom color map
#' color_map <- c(A="#B2DF8A", B="#E31A1C", C="#C41E4E", D="#FEAAAA")
#' names(color_map) <- c("A", "B", "C", "D")
#'
#' plot_tissue(sim, color_map = color_map)
#'
#' # this function returns `TRUE` for cells in the rectangle
#' # [250,350]x[300,350]
#' highlight_function <- function(row) {
#'   (row[["position_x"]] >= 250 && row[["position_x"]] <= 350
#'    && row[["position_y"]] >= 300 && row[["position_y"]] <= 350)
#' }
#'
#' # plot the tissue highlighting the region [250,350]x[300,350]
#' plot_tissue(sim, highlight_function = highlight_function)
#' @export
#'
plot_tissue <- function(simulation, num_of_bins = 100,
                        before_sample = NULL,
                        at_sample = NULL,
                        plot_next_sample_regions = FALSE,
                        plot_sample_region = TRUE,
                        color_map = NULL,
                        list_all_species = FALSE,
                        highlight_function = NULL) {
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))

  sample_info <- NULL

  # Get the cells
  if (is.null(at_sample)) {
    if (is.null(before_sample)) {
      cells <- simulation$get_cells()
    } else {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% dplyr::filter(.data$name == before_sample)

      selected_samples <- samples %>%
        dplyr::filter(.data$time == sample_info$time[[1]])

      first_sample_at_time <- samples %>%
        dplyr::filter(.data$id == min(selected_samples$id))

      cells <- simulation$get_cells(first_sample_at_time$name[[1]])

      from_id <- first_sample_at_time$id[[1]]
    }
  } else {
    if (is.null(before_sample)) {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% dplyr::filter(.data$name == at_sample)

      from_id <- sample_info$id[[1]]

      cells <- simulation$get_cells(at_sample)
    } else {
      stop(paste("The parameters `at_sample` and `before_sample`",
                 "are mutually exclusive."))
    }
  }

  cells <- cells %>% add_species_col()

  if (is.null(color_map)) {
    color_map <- get_species_colors(simulation)
  }

  cells$species <- factor(cells$species,
                          levels = unique(names(color_map)))

  if (is.null(highlight_function)) {
    highlighted_cells <- cells
    masked_cells <- cells[0, ]
  } else {
    highlighted_cells <- cells[apply(cells, 1, highlight_function), ]

    masking_function <- function(row) {
      !highlight_function(row)
    }

    masked_cells <- cells[apply(cells, 1, masking_function), ]
  }

  bin_width_x <- (max(cells$position_x) - min(cells$position_x)) / num_of_bins
  bin_width_y <- (max(cells$position_y) - min(cells$position_y)) / num_of_bins

  pl <- ggplot2::ggplot() +
    ggplot2::geom_hex(
      data = masked_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   fill = .data$species),
      alpha = 0.5,
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = FALSE
    ) +
    ggplot2::geom_hex(
      data = masked_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   fill = "grey70"),
      alpha = 0.5,
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = FALSE
    ) +
    ggplot2::geom_hex(
      data = highlighted_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   fill = .data$species),
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = color_map,
                               drop = !list_all_species) +
    my_theme() +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = "Species") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlim(-10, simulation$get_tissue_size()[1] + 10) +
    ggplot2::ylim(-10, simulation$get_tissue_size()[2] + 10)

  if (!is.null(sample_info)) {
    if (plot_sample_region) {
      pl <- pl + ggplot2::annotate("rect",
                                   xmin = sample_info$xmin[[1]],
                                   xmax = sample_info$xmax[[1]],
                                   ymin = sample_info$ymin[[1]],
                                   ymax = sample_info$ymax[[1]],
                                   fill = NA, color = "black")
    }

    if (plot_next_sample_regions) {
      if (is.null(before_sample)) {
        selected_samples <- samples %>%
          dplyr::filter(.data$time == sample_info$time[[1]],
                        .data$id >= from_id)
      }

      for (row_idx in seq_len(nrow(selected_samples))) {
        sample_info <- selected_samples[row_idx, ]
        pl <- pl + ggplot2::annotate("rect",
                                     xmin = sample_info$xmin[[1]],
                                     xmax = sample_info$xmax[[1]],
                                     ymin = sample_info$ymin[[1]],
                                     ymax = sample_info$ymax[[1]],
                                     fill = NA, color = "black")
      }
    }
  }

  pl
}
