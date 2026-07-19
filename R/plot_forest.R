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

#' Plot a sample forest
#'
#' @description
#' Plot a sample forest. This plot is carried out using
#' `ggraph` and for simplicity of visualisation the forest is plot as a
#' set of trees connected to a generic wildtype cell.
#'
#' @param forest The sample forest to be plot.
#' @param highlight_sample If a sample name, the path from root to the sampled
#'   cells in the sample is highlighted. If `NULL` (default), nothing is
#'   highlighted.
#' @param color_map A named vector representing the simulation species color
#'   map (optional).
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' # use a sample forest example
#' forest <- example("SampleForest")
#'
#' # plot the forest
#' plot_forest(forest)
#'
#' # define a custom color map for the forest species
#' color_map <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99")
#'
#' plot_forest(forest, color_map = color_map)
plot_forest <- function(forest, highlight_sample = NULL, color_map = NULL) {
  if (!inherits(forest, "Rcpp_SampleForest")
        && !inherits(forest, "Rcpp_PhylogeneticForest")) {
    stop("The parameter \"forest\" is not a ProCESS forest.")
  }

  forest_data <- forest$get_nodes()

  if (nrow(forest_data) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest_data %>%
      dplyr::as_tibble() %>%
      dplyr::rename(
        from = .data$ancestor,
        to = .data$cell_id
      ) %>%
      add_species_col() %>%
      dplyr::select(.data$from, .data$to, .data$species,
                    .data$sample, .data$birth_time)

    first_cell <- forest_data %>%
      dplyr::filter(.data$birth_time == 0)

    forest_data <- forest_data %>%
      dplyr::add_row(from = NA, to = NA,
                     species = first_cell[1, ]$species,
                     sample = NA, birth_time = 0) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT", .data$from),
        to = ifelse(is.na(.data$to), "WT", .data$to),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
        highlight = FALSE
      )

    # Highlight is optional
    if (!is.null(highlight_sample)) {
      highlight <- paths_to_sample(forest_data, highlight_sample)
      forest_data$highlight <- forest_data$to %in% highlight$to
    }

    edges <- forest_data %>%
      dplyr::select("from", "to", "highlight")

    # Create a tidygraph object
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)

    graph <- graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(
        forest_data %>%
          dplyr::rename(name = .data$to) %>%
          dplyr::mutate(name = as.character(.data$name)),
        by = "name"
      )

    # Define the layout with the tree rooted and the root on top
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")

    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time

    layout$y <- layout$reversed_btime

    if (is.null(color_map)) {
      color_map <- get_species_colors(forest)
    }

    nsamples <- forest$get_samples_info() %>% nrow()

    labels_every <- max_Y / 10

    point_size <- c(.5, rep(1, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>%
                             dplyr::pull(.data$name))

    # Plot the graph
    graph_plot <- ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = .1,
                             ggplot2::aes(edge_color = ifelse(highlight,
                                                              "indianred3",
                                                              "black")))

    graph_plot +
      ggraph::geom_node_point(ggplot2::aes(color = .data$species,
                                           shape = ifelse(is.na(.data$sample),
                                                          "N/A",
                                                          .data$sample),
                                           size = .data$sample)) +
      ggplot2::scale_shape_manual(values = c(0:nsamples + 1)) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        color = "Species",
        shape = "Sample",
        x = NULL,
        y = "Time"
      ) +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample"),
                      color = ggplot2::guide_legend("Species")) +
      ggplot2::scale_size_manual(
        values = point_size
      ) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y,
                                               labels_every) %>%
          round %>%
          rev,
        breaks = seq(0, max_Y, labels_every) %>% round
      ) +
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
}


paths_to_sample <- function(forest_data, sample) {
  to <- forest_data %>% dplyr::filter(sample == !!sample) %>%
    dplyr::pull(to)

  queue <- NULL

  while (length(to) > 0) {
    # New element to inspect
    to_head <- to[1]
    to <- to[-1]

    # Forward star
    to_add <- forest_data %>%  dplyr::filter(to == to_head)

    if (to_add$from != to_head) {
      queue <- dplyr::bind_rows(queue, to_add) %>% dplyr::distinct()

      # Recursion
      to <- c(to, to_add %>% dplyr::pull(.data$from)) %>% unique()
    }
  }

  queue
}

