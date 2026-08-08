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

sample_shape <- function(nodes) {
  nodes %>%
    dplyr::mutate(shape_label = ifelse(is.na(.data$sample), "N/A",
                                       .data$sample)) %>%
    dplyr::pull(.data$shape_label)
}

#' Plot a sample forest
#'
#' @description
#' Plot a sample forest. This plot is carried out using
#' `ggraph` and for simplicity of visualisation the forest is plot as a
#' set of trees connected to a generic wild-type cell.
#'
#' @param forest The sample forest to be plot.
#' @param highlight_sample If a sample name, the path from root to the sampled
#'   cells in the sample is highlighted. If `NULL` (default), nothing is
#'   highlighted.
#' @param color_map A named vector representing the simulation species color
#'   map (optional when `color_label_function` is `NULL`; mandatory,
#'   otherwise).
#' @param alpha_function A function whose input is the data frame
#'   returned by <code>[SampleForest$get_cells()]</code> or
#'   <code>[PhylogeneticForest$get_cells()]</code> and returns
#'   a real vector whose values are in the interval \eqn{[0,1]} and whose
#'   length is the number of rows in the input data frame. Each value in the
#'   output is used as alpha level of the corresponding cell. When the
#'   parameter is set to `NULL`, all nodes have alpha level `1`
#'   (default: `NULL`).
#' @param shape_label_function A function whose input is the data frame
#'   returned by <code>[SampleForest$get_cells()]</code> or
#'   <code>[PhylogeneticForest$get_cells()]</code> and returns
#'   a string vector whose length is the number of rows in the input data
#'   frame. Each value in the output is the label of the associated node
#'   and determine the node shape. When the parameter is set to `NULL`,
#'   all nodes have the same shape. By default, the node shapes correspond
#'   to the sample which collected the corresponding cell.
#' @param color_label_function A function whose input is the data frame
#'   returned by <code>[SampleForest$get_cells()]</code> or
#'   <code>[PhylogeneticForest$get_cells()]</code> and returns
#'   a string vector whose length is the number of rows in the input data
#'   frame. Each value in the output is the label of the associated node
#'   and determine the node color. When the parameter is set to `NULL`,
#'   the nodes are colored according their species (default: `NULL`).
#' @return A `ggraph` tree plot.
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
#'
#' # define an alpha function assigning to each node an alpha
#' # value according to the corresponding cell's birth time
#' alpha_f <- function(nodes) {
#'   library(dplyr)
#'
#'   T <- max(nodes$birth_time)
#'   nodes %>%
#'     mutate(alpha = (T-.data$birth_time)/T) %>%
#'     pull(alpha)
#' }
#'
#' plot_forest(forest, alpha_function = alpha_f)
#'
#' # plot the forest avoiding to denote samples by using node shapes
#' plot_forest(forest, shape_label_function = NULL)
#'
#' # a shape function that labels each node as the sample that
#' # collected the corresponding cells
#' shape_label_function <- function(nodes) {
#'   library(dplyr)
#'
#'   nodes %>%
#'     mutate(label = case_when(birth_time <= 300 ~ "Old",
#'                              TRUE ~ "Young")) %>%
#'     pull(label)
#' }
#'
#' # plot the forest using the node shapes to denote cell birth time
#' plot_forest(forest, shape_label_function = shape_label_function)
#'
#' color_label_function <- function(nodes) {
#'   library(dplyr)
#'
#'   nodes %>%
#'     mutate(age = case_when(birth_time <= 300 ~ "Old",
#'                            TRUE ~ "Young")) %>%
#'     mutate(label = paste(mutant, age)) %>%
#'     pull(label)
#' }
#'
#' # get the plot labels (i.e., paste(mutants, "Old") +
#' # paste(mutants, "Young"))
#' mutants <- unique(forest$get_nodes() %>% dplyr::pull(mutant))
#' labels <- c(paste(mutants, "Old"), paste(mutants, "Young"))
#'
#' # create a color map for the labels
#' color_map <- RColorBrewer::brewer.pal(n = length(labels), name = "Set1")
#' names(color_map) <- labels
#'
#' # plot the tissue labelling the cells according to
#' # `label_function`. The parameter `color_map` is mandatory.
#' plot_forest(forest, color_label_function = color_label_function,
#'             color_map = color_map, shape_label_function = NULL)
#' @export
#'
plot_forest <- function(forest, highlight_sample = NULL, color_map = NULL,
                        alpha_function = NULL,
                        shape_label_function = sample_shape,
                        color_label_function = NULL) {
  if (!inherits(forest, "Rcpp_SampleForest")
      && !inherits(forest, "Rcpp_PhylogeneticForest")) {
    stop("The parameter \"forest\" is not a ProCESS forest.")
  }

  forest_data <- forest$get_nodes()

  if (is.null(color_label_function)) {
    forest_data <- forest_data %>% add_species_col("color_label")

    if (is.null(color_map)) {
      color_map <- get_species_colors(forest)
    }
  } else {
    if (is.null(color_map)) {
      stop("\"color_map\" is mandatory when ",
           "\"color_label_function\" is specified.")
    }

    forest_data[["color_label"]] <- color_label_function(forest_data)
  }

  forest_data$color_label <- factor(forest_data$color_label,
                                    levels = names(color_map))

  if (!is.null(shape_label_function)) {
    forest_data[["shape_label"]] <- shape_label_function(forest_data)
  } else {
    forest_data[["shape_label"]] <- NA
  }

  if (is.null(alpha_function)) {
    forest_data <- forest_data %>%
      dplyr::mutate(alpha_level = 1)
  } else {
    forest_data <- forest_data %>%
      dplyr::mutate(alpha_level = alpha_function(.))
  }

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
      dplyr::select(.data$from, .data$to, .data$sample,
                    .data$shape_label, .data$color_label,
                    .data$birth_time, .data$alpha_level)

    first_cell <- forest_data %>%
      dplyr::filter(.data$birth_time == 0)

    forest_data <- forest_data %>%
      dplyr::add_row(from = NA, to = NA,
                     color_label = first_cell[1, ]$color_label,
                     shape_label = first_cell[1, ]$shape_label,
                     birth_time = 0) %>%
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

    graph <- graph %>%
      tidygraph::activate("edges") %>%
      dplyr::mutate(edge_alpha = tidygraph::.N()$alpha_level[to])

    # Define the layout with the tree rooted and the root on top
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")

    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time

    layout$y <- layout$reversed_btime

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
                                                              "black"),
                                          alpha = .data$edge_alpha))

    if (is.null(shape_label_function)) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(color = .data$color_label,
                                             alpha = .data$alpha_level,
                                             size = .data$sample))
    } else {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(color = .data$color_label,
                                             shape = .data$shape_label,
                                             alpha = .data$alpha_level,
                                             size = .data$sample)) +
        ggplot2::scale_shape_manual(values = c(0:nsamples + 1))
    }

    graph_plot +
      ggplot2::scale_alpha_identity() +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        color = NULL,
        shape = NULL,
        x = NULL,
        y = "Time"
      ) +
      ggplot2::guides(size = "none",
                      shape = NULL,
                      color = NULL) +
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
