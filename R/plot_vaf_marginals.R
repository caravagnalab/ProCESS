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

# add labels for driver mutations
add_marginals_driver_mutation_labels <- function(plot, data, driver_mutations) {

  drivers <- data %>% dplyr::filter(.data$classes == "driver")

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
        fill = causes
      ),
      color = 'black',
      size = 2,
      min.segment.length = 0,
      segment.color = "grey50", # Color of the connecting segment
      segment.linetype = "dashed", # Line type of the segment
      nudge_x = 0.1,
      nudge_y = 0.1,
      show.legend = FALSE,
      segment.curvature = 0,
      segment.ncp = 1,
      segment.square = TRUE,
      segment.inflect = TRUE,
      direction = "both"
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
#'   `function(x) x %>% dplyr::filter(classes != "germinal")`).
#' @param driver_mutations The data frame of the driver mutations as
#'   returned by `PhylogeneticForest$get_driver_mutations()`.
#'   This parameter can be avoided when `seq_result` is the result
#'   of the function `simulate_seq()` (default: `NULL`).
#' @param driver_mutation_labels A Boolean value to enable/disable
#'   driver mutation labels in the returned plot (default: TRUE).
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
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.1,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#' sim$sample_cells("SampleA", bbox$lower_corner, bbox$upper_corner)
#'
#' # adding second mutant
#' sim$add_mutant(name = "B",
#'                growth_rates = 0.3,
#'                death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # sampling tissue again
#' bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
#' sim$sample_cells("SampleB", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_sample_forest()
#'
#' # placing mutations
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8),
#'                     drivers = list(SNV("22", 16510210, "C", "T", allele = 1),
#'                                    "DGCR8 P26L"))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9),
#'                     drivers = list("DGCR8 A18V"))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' # simulating sequencing without the normal sample
#' seq_results <- simulate_seq(phylo_forest, coverage = 10, write_SAM = F,
#'                             with_normal_sample = FALSE, quiet = TRUE)
#'
#' # plotting the VAF marginals without germinal mutations
#' plot_VAF_marginals(seq_results)
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!classes %in% list("germinal",
#'                                             "pre-neoplastic"))
#' }
#'
#' # plotting the VAF marginals without germinal and pre-neoplastic
#' plot_VAF_marginals(seq_results, mutation_filter=filter_data)
#'
#' # plotting the VAF marginals filtering out mutations having
#' # the VAF below 0.2 in both the samples
#' plot_VAF_marginals(seq_results, cuts = c(0.2, 1))
#'
#' # plotting the VAF marginals with labels
#' plot_VAF_marginals(seq_results, labels = seq_results$mutations["causes"])
#'
#' # avoid the driver mutation labels
#' plot_VAF_marginals(seq_results, labels = seq_results$mutations["causes"],
#'                    driver_mutation_labels = FALSE)
#'
#' # the same plots can be drawn by using the mutations data frame
#' # in place of the `simulate_seq()` output
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>%
#'    dplyr::filter(classes!="germinal")
#'
#' # plotting the VAF histogram filtering out VAFs below 0.02
#' plot_VAF_marginals(f_seq, cuts = c(0.2, 1))
#'
#' # plotting the VAF histogram with labels
#' plot_VAF_marginals(f_seq, labels = f_seq["causes"])
#'
#' # use the driver codes in the driver mutation labels
#' plot_VAF_marginals(f_seq, labels = f_seq["causes"],
#'                    driver_mutations = phylo_forest$get_driver_mutations())
#'
#' # avoid the driver mutation labels
#' plot_VAF_marginals(f_seq, labels = f_seq["causes"],
#'                    driver_mutation_labels = FALSE)
#'
#' # deleting the mutation engine directory
#' unlink('demo', recursive = T)
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
                                             "classes", "causes")) %>%
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
                                               col = .data$labels.x))
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
