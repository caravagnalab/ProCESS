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

#' Plot a Variant Allele Frequency (VAF) histogram
#'
#' This function generates a histogram showing the distribution of Variant
#' Allele Frequency (VAF) across samples and chromosomes.
#'
#' @param seq_result A data frame containing sequencing results.
#' @param chromosomes A character vector specifying the chromosomes to
#'   include in the plot (default: all the chromosomes in `seq_res`).
#' @param samples A character vector specifying the sample names to include
#'   in the plot. When set to `NULL`, the function includes all samples
#'   except the "normal_sample" (default: `NULL`).
#' @param colour_by A character indicating whether to color the histogram
#'   bars by "causes" or "classes" (default: "causes").
#' @param binwidth The width of the plot bins. When set to `NULL`, the
#'   function computes the most convenient bin width according to the
#'   maximum coverage reported in the data frame (default: `NULL`).
#' @param cuts A numeric vector specifying the range of VAF values to
#'   include in the plot (default: `c(0, 1)`).
#' @return A ggplot2 object showing the VAF histogram.
#' @seealso `plot_VAF_marginals()`, `plot_VAF()`
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
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' # simulating sequencing without the normal sample
#' seq_results <- simulate_seq(phylo_forest, coverage = 10, write_SAM = F,
#'                             with_normal_sample = FALSE)
#'
#' library(dplyr)
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>% dplyr::filter(classes!="germinal")
#'
#' # plotting the VAF histogram
#' plot_VAF_histogram(f_seq, cuts = c(0.02, 1))
#'
#' # plotting the VAF histogram with labels
#' plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1))
#'
#' # deleting the mutation engine directory
#' unlink('demo', recursive = T)
plot_VAF_histogram <- function(
    seq_result,
    chromosomes = NULL,
    samples = NULL,
    labels = NULL,
    binwidth = NULL,
    cuts = c(0, 1)
) {
  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_result) && ("mutations" %in% names(seq_result))) {

    # extract the field
    seq_res <- seq_result["mutations"]
  } else {
    seq_res <- seq_result
  }

  data <- seq_to_long(seq_res)

  if (!is.null(labels)) {
    if (!is(labels, "data.frame")) {
      stop("The parameter \"labels\" must be a data frame when non-NULL.")
    }

    if (length(labels) != 1) {
      stop(paste0("The parameters \"labels\" must be a data frame ",
                  "with one column when non-NULL."))
    }

    if (nrow(labels) != nrow(seq_res)) {
      stop(paste0("The parameters \"seq_res\" and \"labels\"",
                  " must have the same number of rows."))
    }

    data["labels"] <- labels
    label_name <- names(labels)
  }

  chromosomes <- validate_chromosomes(seq_res, chromosomes)

  data <- data %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample_name))) {
      stop("Invalid sample name in samples parameter")
    }
    data <- data %>% dplyr::filter(sample_name %in% samples)
  }

  if (!is.null(labels)) {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF,
                                             fill = labels)) +
      ggplot2::labs(col = data$labels, fill = label_name)
  } else {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF))
  }

  if (is.null(binwidth)) {
    binwidth <- 1 / max(data$DP)
  }

  plot +
    ggplot2::geom_histogram(binwidth = binwidth, alpha = 0.5) +
    ggplot2::facet_grid(sample_name ~ chr, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    ggplot2::theme(legend.position = "bottom")
}