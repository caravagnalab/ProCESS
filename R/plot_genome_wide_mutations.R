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

#' Plot the genome-wide Depth Ratio (DR)
#'
#' @description This function plots the genome-wide Depth Ratio (DR)
#'   of a specific sample.
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
#' @param cuts A numeric vector specifying the range of VAF values to
#'   include in the plot (default: `c(0, 1)`).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the DR distribution across the genome.
#' @seealso `plot_VAF()`, `plot_BAF()`
#' @examples
#' # use a sequencing result example
#' seq_results <- example("Sequencing results")
#'
#' # plotting the depth ratio over all samples
#' plot_DR(seq_results)
#'
#' # plotting the depth ratio for sample S_2_2 only
#' plot_DR(seq_results, samples = "S_2_2")
#'
#' # plotting the depth ratio for S_2_2 with labels
#' plot_DR(seq_results, samples = "S_2_2",
#'         labels = seq_results$mutations["nature"])
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!nature %in% list("germinal",
#'                                            "pre-neoplastic"))
#' }
#'
#' # plotting the depth ratio without germinal and pre-neoplastic
#' plot_DR(seq_results, mutation_filter = filter_data)
#'
#' # the same plots can be drawn by using the mutations data frame
#' # in place of the `simulate_seq()` output
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>%
#'    dplyr::filter(nature != "germinal")
#'
#' # plotting the depth ratio
#' plot_DR(f_seq)
#'
#' @export
#'
plot_DR <- function(
    seq_result,
    chromosomes = NULL,
    samples = NULL,
    labels = NULL,
    mutation_filter = filter_germinal,
    cuts = c(0, 1),
    N = 5000) {

  setup_data <- setup_data_for_plotting(seq_result, chromosomes,
                                        samples, labels, mutation_filter,
                                        NULL)

  if (nrow(setup_data$normal) == 0) {
    stop(paste("The parameter \"seq_result\" does not contain",
               "the mandatory normal sample \"normal sample\"."))
  }

  data <- setup_data$data %>%
    dplyr::filter(.data$sample != "normal.sample")

  N <- min(N, nrow(data))

  data <- data %>%
    dplyr::left_join(setup_data$normal, suffix = c(".tumour", ".normal"),
                     by = c("chr", "from", "ref", "alt")) %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>%
    dplyr::sample_n(N) %>%
    dplyr::arrange(.data$chr, .data$from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n())

  chr_limits <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::filter(.data$abs_pos == min(.data$abs_pos)) %>%
    dplyr::pull(.data$abs_pos)
  chr_limits <- c(chr_limits, max(data$abs_pos))

  chr_means <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(mean_pos = (max(.data$abs_pos) +
                                   min(.data$abs_pos)) / 2) %>%
    dplyr::pull(.data$mean_pos)

  if (!is.null(labels)) {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$DR,
                                             color = .data$labels.tumour)) +
      ggplot2::labs(color = setup_data$label_name)
  } else {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$DR))
  }

  plot +
    ggplot2::facet_wrap(~sample.tumour,  ncol = 1,
                        strip.position = "right") +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(x = NULL, y = "DR") +
    ggplot2::lims(y = c(0, NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means,
                                labels = unique(data$chr))
}


#' Plot the genome-wide B-Allele Frequency (BAF)
#'
#' @description This function plots the genome-wide B-Allele Frequency (BAF)
#'   of a specific sample.
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
#'   include in the plot (default: `c(0, 1)`).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the BAF distribution across the genome.
#' @seealso `plot_VAF()`, `plot_DR()`
#' @examples
#' # use a sequencing result example
#' seq_results <- example("Sequencing results")
#'
#' # plotting the BAF over all samples
#' plot_BAF(seq_results)
#'
#' # plotting the BAF for sample S_2_2 only
#' plot_BAF(seq_results, samples = c("S_2_2"))
#'
#' # plotting the BAF for S_2_2 with labels
#' plot_BAF(seq_results, samples = "S_2_2",
#'          labels = seq_results$mutations["nature"])
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!nature %in% list("germinal",
#'                                            "pre-neoplastic"))
#' }
#'
#' # plotting the BAF without germinal and pre-neoplastic
#' plot_BAF(seq_results, mutation_filter = filter_data)
#'
#' # the same plots can be drawn by using the mutations data frame
#' # in place of the `simulate_seq()` output
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>%
#'    dplyr::filter(nature != "germinal")
#'
#' # plotting the BAF
#' plot_BAF(f_seq)
#'
#' @export
#'
plot_BAF <- function(
    seq_result,
    chromosomes = NULL,
    samples = NULL,
    labels = NULL,
    mutation_filter = filter_germinal,
    driver_mutations = NULL,
    driver_mutation_labels = TRUE,
    cuts = c(0, 1),
    N = 5000) {
  setup_data <- setup_VAF_data_for_plotting(seq_result, chromosomes,
                                            samples, labels, mutation_filter,
                                            driver_mutations, cuts)

  tumour_data <- setup_data$data %>%
    dplyr::filter(.data$sample != "normal.sample")

  tumour_data$BAF <- pmin(tumour_data$VAF, 1 - tumour_data$VAF)

  tumour_data_no_drivers <- tumour_data %>%
    dplyr::filter(.data$nature != "driver")

  N <- min(N, nrow(tumour_data_no_drivers))

  data <- tumour_data_no_drivers %>%
    dplyr::sample_n(N)

  driver_data <- tumour_data %>%
    dplyr::filter(.data$nature == "driver")

  data <- dplyr::bind_rows(data, driver_data) %>%
    dplyr::arrange(.data$chr, .data$from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n())

  chr_limits <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::filter(.data$abs_pos == min(.data$abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(data$abs_pos))

  chr_means <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(mean_pos = (max(.data$abs_pos) +
                                   min(.data$abs_pos)) / 2) %>%
    dplyr::pull(.data$mean_pos)

  if (!is.null(labels)) {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$BAF,
                                             color = labels)) +
      ggplot2::labs(color = setup_data$label_name)
  } else {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$BAF))
  }

  plot <- plot +
    ggplot2::facet_wrap(~sample,  ncol = 1, strip.position = "right") +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(x = NULL, y = "BAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(data$chr))

  if (driver_mutation_labels) {
    plot <- add_driver_mutation_labels(plot, data, setup_data$driver_mutations,
                                       "abs_pos", "BAF")
  }

  plot
}


#' Plot the genome-wide Variant Allele Frequency (VAF)
#'
#' @description This function plots the genome-wide Variant Allele
#'   Frequency (VAF) of a specific sample.
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
#'   include in the plot (default: `c(0, 1)`).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the BAF distribution across the genome.
#' @seealso `plot_BAF()`, `plot_DR()`, `plot_VAF_histogram()`
#' @examples
#' # use a sequencing result example
#' seq_results <- example("Sequencing results")
#'
#' # plotting the VAF over all samples
#' plot_VAF(seq_results)
#'
#' # plotting the VAF for sample S_2_2 only
#' plot_VAF(seq_results, samples = c("S_2_2"))
#'
#' # plotting the VAF for S_2_2 with labels
#' plot_VAF(seq_results, samples = "S_2_2",
#'          labels = seq_results$mutations["nature"])
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!nature %in% list("germinal",
#'                                            "pre-neoplastic"))
#' }
#'
#' # plotting the VAF without germinal and pre-neoplastic
#' plot_VAF(seq_results, mutation_filter = filter_data)
#'
#' # the same plots can be drawn by using the mutations data frame
#' # in place of the `simulate_seq()` output
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>%
#'    dplyr::filter(nature != "germinal")
#'
#' # plotting the VAF
#' plot_VAF(f_seq)
#'
#' @export
#'
plot_VAF <- function(
    seq_result,
    chromosomes = NULL,
    samples = NULL,
    labels = NULL,
    mutation_filter = filter_germinal,
    driver_mutations = NULL,
    driver_mutation_labels = TRUE,
    cuts = c(0, 1),
    N = 5000) {
  setup_data <- setup_VAF_data_for_plotting(seq_result, chromosomes,
                                            samples, labels, mutation_filter,
                                            driver_mutations, cuts)

  tumour_data <- setup_data$data %>%
    dplyr::filter(.data$sample != "normal.sample")

  tumour_data_no_drivers <- tumour_data %>%
    dplyr::filter(.data$nature != "driver")

  N <- min(N, nrow(tumour_data_no_drivers))

  data <- tumour_data_no_drivers %>%
    dplyr::sample_n(N)

  driver_data <- tumour_data %>%
    dplyr::filter(.data$nature == "driver")

  data <- dplyr::bind_rows(data, driver_data) %>%
    dplyr::arrange(.data$chr, .data$from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n())

  chr_limits <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::filter(.data$abs_pos == min(.data$abs_pos)) %>%
    dplyr::pull(.data$abs_pos)
  chr_limits <- c(chr_limits, max(data$abs_pos))

  chr_means <- data %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(mean_pos = (max(.data$abs_pos) +
                                   min(.data$abs_pos)) / 2) %>%
    dplyr::pull(.data$mean_pos)

  if (!is.null(labels)) {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$VAF,
                                             color = labels)) +
      ggplot2::labs(color = setup_data$label_name)
  } else {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$abs_pos,
                                             y = .data$VAF))
  }

  plot <- plot +
    ggplot2::facet_wrap(~sample,  ncol = 1, strip.position = "right") +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(x = NULL, y = "VAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(data$chr))

  if (driver_mutation_labels) {
    plot <- add_driver_mutation_labels(plot, data, setup_data$driver_mutations,
                                       "abs_pos", "VAF")
  }

  plot
}
