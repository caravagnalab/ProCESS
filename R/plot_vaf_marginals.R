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

# add labels for driver mutations
add_marginals_driver_mutation_labels <- function(plot, data, driver_mutations) {

  drivers <- data %>% dplyr::filter(.data$nature == "driver")

  if (!"mutant" %in% colnames(drivers)) {
    drivers <- drivers %>%
      dplyr::mutate(mutant = .data$cause)
  }

  if (!is.null(driver_mutations)) {

    drivers <- dplyr::left_join(drivers, driver_mutations,
                                by = c("chr" = "chr",
                                       "from" = "start",
                                       "ref" = "ref",
                                       "alt" = "alt"))
  } else {
    drivers$code <- rep(NA, times = nrow(drivers))
  }

  for (i in seq_len(nrow(drivers))) {
    driver <- drivers[i,]
    if (is.na(driver[["code"]])) {
      driver[["code"]] <- mutation_string(driver)
    }
    drivers[i,] <- driver
  }

  plot <- plot +
    ggrepel::geom_label_repel(
      data = drivers,
      ggplot2::aes(
        x = VAF.x,
        y = VAF.y,
        label = code,
        fill = mutant
      ),
      color = 'black',
      size = 2,
      min.segment.length = 0,
      segment.color = "grey50", # Color of the connecting segment
      segment.linetype = "dashed", # Line type of the segment
      key_glyph = "rect",
      nudge_x = 0.1,
      nudge_y = 0.1,
      #show.legend = FALSE,
      segment.curvature = 0,
      segment.ncp = 1,
      segment.square = TRUE,
      segment.inflect = TRUE,
      direction = "both",
      box.padding = 0.5
    )

  return(plot)
}

#' Plot Marginals of Variant Allele Frequency (VAF)
#'
#' This function generates scatter plots showing the marginal distributions of
#' Variant Allele Frequency (VAF) for pairs of samples on a specific chromosome.
#'
#' @param seq_result Either the output of `simulate_seq()` or a data
#'   frame containing sequencing results.
#' @param chromosomes A character vector specifying the chromosomes to
#'   include in the plot (default: all the chromosomes in `seq_res`).
#' @param samples A character vector specifying the sample names to include
#'   in the plot. When set to `NULL`, the function includes all samples
#'   except the "normal sample" (default: `NULL`).
#' @param labels A data frame column labelling each mutation. Each label
#'   is associated to a different colour in the plot (default: `NULL`).
#' @param mutation_filter A function filtering mutations from the input
#'   data (default: a function filtering out "germinal" mutations, e.g.,
#'   `function(x) x %>% dplyr::filter(nature != "germinal")`).
#' @param driver_mutations The data frame of the driver mutations as
#'   returned by `PhylogeneticForest$get_driver_mutations()`.
#'   This parameter can be avoided when `seq_result` is the result
#'   of the function `simulate_seq()` (default: `NULL`).
#' @param driver_mutation_labels A Boolean value to enable/disable
#'   driver mutation labels in the returned plot (default: `TRUE`).
#' @param cuts A numeric vector specifying the range of VAF values to
#'   include in the plot. A mutation is represented in the plot only if
#'   its VAF lays in this interval in at least one of the samples
#'   (default: `c(0, 1)`).
#' @return A list of ggplot2 objects showing scatter plots of VAF marginals
#'   for pairs of samples.
#' @seealso `plot_VAF_histogram()`, `plot_VAF()`
#' @export
#'
#' @examples
#' # use a sequencing result example
#' seq_results <- example("Sequencing results")
#'
#' # plotting the VAF marginals without germinal mutations
#' plot_VAF_marginals(seq_results)
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!nature %in% list("germinal",
#'                                            "pre-neoplastic"))
#' }
#'
#' # plotting the VAF marginals without germinal and pre-neoplastic
#' plot_VAF_marginals(seq_results, mutation_filter = filter_data)
#'
#' # plotting the VAF marginals filtering out mutations having
#' # the VAF below 0.2 in both the samples
#' plot_VAF_marginals(seq_results, cuts = c(0.2, 1))
#'
#' # plotting the VAF marginals with labels
#' plot_VAF_marginals(seq_results, labels = seq_results$mutations["cause"])
#'
#' # avoid the driver mutation labels
#' plot_VAF_marginals(seq_results, labels = seq_results$mutations["cause"],
#'                    driver_mutation_labels = FALSE)
#'
#' # the same plots can be drawn by using the mutations data frame
#' # in place of the `simulate_seq()` output
#'
#' # get driver mutations
#' driver_mutations <- seq_results$parameters$driver_mutations
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>%
#'    dplyr::filter(nature!="germinal")
#'
#' # plotting the VAF histogram filtering out VAFs below 0.02
#' plot_VAF_marginals(f_seq, cuts = c(0.2, 1))
#'
#' # plotting the VAF histogram with labels
#' plot_VAF_marginals(f_seq, labels = f_seq["cause"])
#'
#' # use the driver codes in the driver mutation labels
#' plot_VAF_marginals(f_seq, labels = f_seq["cause"],
#'                    driver_mutations = driver_mutations)
#'
#' # avoid the driver mutation labels
#' plot_VAF_marginals(f_seq, labels = f_seq["cause"],
#'                    driver_mutation_labels = FALSE)
#'
plot_VAF_marginals <- function(
    seq_result,
    chromosomes = NULL,
    samples = NULL,
    labels = NULL,
    mutation_filter = filter_germinal,
    driver_mutations = NULL,
    driver_mutation_labels = TRUE,
    cuts = c(0, 1))
{
  setup_data <- setup_VAF_data_for_plotting(seq_result, chromosomes, samples,
                                            labels, mutation_filter,
                                            driver_mutations, cuts)

  data <- setup_data$data
  driver_mutations <- setup_data$driver_mutations

  if (length(unique(data$sample)) < 2) {
    stop("At least two samples are required.")
  }

  combinations <- utils::combn(unique(data$sample), m = 2)

  lapply(seq_len(ncol(combinations)), function(i) {
    couple <- combinations[, i]
    d1 <- data %>% dplyr::filter(.data$sample == couple[1])
    d2 <- data %>% dplyr::filter(.data$sample == couple[2])

    djoin <- dplyr::full_join(d1, d2, by = c("chr", "from", "ref", "alt",
                                             "nature", "cause")) %>%
      dplyr::mutate(VAF.x = ifelse(is.na(.data$VAF.x), 0, .data$VAF.x)) %>%
      dplyr::mutate(VAF.y = ifelse(is.na(.data$VAF.y), 0, .data$VAF.y))

    if (!is.null(labels)) {
      djoin <- djoin %>%
        dplyr::mutate(labels.x = ifelse(is.na(.data$labels.x),
                                        .data$labels.y, .data$labels.x))
    }

    if (!is.null(labels)) {
      plot <- djoin %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = .data$VAF.x,
                                               y = .data$VAF.y,
                                               col = .data$labels.x)) +
        ggplot2::labs(col = djoin$labels.x)
    } else {
      plot <- djoin %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = .data$VAF.x,
                                               y = .data$VAF.y))
    }

    plot <- plot +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::xlim(c(-0.01, 1.01)) +
      ggplot2::ylim(c(-0.01, 1.01)) +
      ggplot2::labs(x = couple[1], y = couple[2]) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    if (driver_mutation_labels) {
      plot <- add_marginals_driver_mutation_labels(plot, djoin,
                                                   driver_mutations)
    }

    plot
  })
}
