## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2026 - Alberto Casagrande <alberto.casagrande@uniud.it>
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

#' Annotate a plot of cell divisions
#'
#' @description
#' It annotates a plot of cell divisions where branches containing relevant
#'   biological events are colored
#'
#' @param forest The original forest object has been derived.
#' @param labels A data frame annotating the sticks (it can be the output of
#'   `get_relevant_branches()`). If no labels are provided the output of
#'   `get_relevant_branches()` is used by default.
#' @param cls A custom list of colors for any stick. If NULL a default palette
#'   is chosen.
#'
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' sim <- TissueSimulation()
#'
#' sim$add_mutant("A", c(duplication = 1))
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_size("A",1e4)
#'
#' sim$add_mutant("B", c(duplication = 3.5))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B",1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100, "B" = 100), 50, 50)
#' sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
#' forest <- sim$get_sample_forest()
#'
#' labels <- get_relevant_branches(forest)
#' plot_sticks(forest, labels)
#' @examples
#' # build a simulation with epigenetic states
#' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
#'
#' # add mutant "A" and set its species rates
#' sim$add_mutant("A",
#'                list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
#'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
#'
#' # add mutant "B" and set its species rates
#' sim$add_mutant("B",
#'                list(E1 = list(duplication = 0.3, death = 0.1, E2 = 0.02),
#'                     E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))
#'
#' # schedule a mutation from "A" to "B"
#' sim$schedule_mutation("A", "B", 20)
#'
#' # place an "A[E1]" cell in the tissue
#' sim$place_cell("A[E1]", 500, 500)
#'
#' # run the simulation up to time 70
#' sim$run_up_to_time(70)
#'
#' # sample the tissue
#' bbox <- sim$search_sample(c("A" = 100, "B" = 100), 50, 50)
#' sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
#'
#' # build the sample forest
#' forest <- sim$get_sample_forest()
#'
#' plot_sticks(forest)

plot_sticks <- function(forest, labels = NULL, cls = NULL) {
  if (!inherits(forest, "Rcpp_SampleForest") &&
      !inherits(forest, "Rcpp_PhylogeneticForest")) {
    stop("The parameter must be a forest.")
  }

  forest_data <- forest$get_nodes() %>%
    add_species_col() %>%
    dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
    dplyr::select(.data$from, .data$to, .data$species,
                  .data$sample, .data$birth_time)

  if (nrow(forest_data) == 0) {
    stop("The forest does not contain any node")
  } else {
    if (is.null(labels)) {
      labels <- get_relevant_branches(forest)
    }

    labels <- labels %>%
      dplyr::select(.data$cell_id, .data$label) %>%
      dplyr::add_row(cell_id = NA, label = "Truncal")

    forest_data <- forest_data %>%
      dplyr::add_row(from = NA, to = NA, species = "Wild-type",
                     sample = NA, birth_time = 0) %>%
      dplyr::full_join(labels, by = dplyr::join_by(x$to == y$cell_id)) %>%
      dplyr::mutate(
        label = ifelse(is.na(.data$label), "Subclonal", .data$label),
        from = ifelse(is.na(.data$from), "WT", .data$from),
        to = ifelse(is.na(.data$to), "WT", as.character(.data$to)),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
        highlight = FALSE
      )

    edges <- forest_data %>% dplyr::select("from", "to", "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>% tidygraph::activate("nodes") %>%
      dplyr::left_join(forest_data %>% dplyr::rename(name = .data$to) %>%
                         dplyr::mutate(name = as.character(.data$name)),
                       by = "name")
    layout <- ggraph::create_layout(graph, layout = "tree",
                                    root = "WT")
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time
    layout$y <- layout$reversed_btime

    if (is.null(cls)) {
      cls <- get_colors_for(forest_data %>%
                              dplyr::arrange(dplyr::desc(.data$to)) %>%
                              dplyr::pull(.data$label) %>% unique())
      cls["Subclonal"] <- "gainsboro"
    }
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y/10
    not_subclonal <- forest_data$label[forest_data$label != "Subclonal"] %>%
      unique()
    point_size <- c(0, rep(1, length(not_subclonal)))
    names(point_size) <- c("Subclonal", not_subclonal)

    ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.5,
                             ggplot2::aes(edge_color = ifelse(highlight,
                                                              "indianred3",
                                                              "gainsboro"))) +

      ggraph::geom_node_point(ggplot2::aes(
        color = .data$label,
        size = .data$label
      )) +
      ggplot2::scale_color_manual(values = cls, na.translate = FALSE)  +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = NULL, y = "Birth time") +
      ggplot2::guides(size = "none",
        shape = ggplot2::guide_legend("Sample"),
        fill = ggplot2::guide_legend("Species")
      ) +
      ggplot2::scale_size_manual(values = point_size) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>%
                                    round %>% rev,
                                  breaks = seq(0, max_Y, labels_every) %>%
                                    round) +
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
}
