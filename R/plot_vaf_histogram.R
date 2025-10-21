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

# compute the bin width for a set of data according to the
# Freedman-Diaconis rule
freedman_diaconis_binwidth <- function(data) {
  q1 <- quantile(data, 0.25)[[1]]
  q3 <- quantile(data, 0.75)[[1]]

  IQR <- q3 - q1

  return(2 * (IQR / ((length(data))^(1 / 3))))
}

# compute a string representation for a mutation row
mutation_string <- function(mutation) {
  paste0(mutation[["chr"]], ":", mutation[["from"]], "[",
         mutation[["ref"]], ">", mutation[["alt"]], "]")
}

# add labels for driver mutations
add_driver_mutation_labels <- function(
    plot,
    data,
    driver_mutations,
    x_value,
    y_value) {

  ggb <- ggplot2::ggplot_build(plot)

  facet_p <- ggb$layout$layout %>% dplyr::mutate(y_max = NA)
  for (i in seq_len(nrow(facet_p))) {
    facet_p$y_max[i] <- ggb$layout$panel_params[[i]]$y.range[2]
  }

  drivers <- data %>% dplyr::filter(.data$classes == "driver") %>%
    dplyr::left_join(facet_p %>% dplyr::select(.data$sample,
                                               .data$y_max,
                                               .data$PANEL),
                     by = c("sample"))

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

  panels <- unique(drivers$PANEL)

  get_value_expr <- function(x) {
    if (is.character(x)) {
      rlang::sym(x)
    } else {
      x
    }
  }

  x_expr <- get_value_expr(x_value)
  y_expr <- get_value_expr(y_value)

  for (panel in panels) {
    drivers_in_panel <- drivers %>% dplyr::filter(PANEL == panel)

    y_max_count <- drivers_in_panel[["y_max"]][1]

    plot <- plot +
      ggrepel::geom_label_repel(
        data = drivers_in_panel,
        ggplot2::aes(
          x = !!x_expr,
          y = !!y_expr,
          label = code,
          fill = causes
        ),
        color = 'black',
        ylim = c(y_max_count * 0.5,
                 y_max_count * 0.9),
        size = 2,
        min.segment.length = 0,
        segment.color = "grey50", # Color of the connecting segment
        segment.linetype = "dashed", # Line type of the segment
        show.legend = FALSE,
        segment.curvature = 0,
        segment.ncp = 1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        direction = "both"
      )
  }

  return(plot)
}

# filter germinal mutations from data
filter_germinal <- function(data) {
  data %>% dplyr::filter(.data$classes != "germinal")
}

# setup input data for plotting
setup_data_for_plotting <- function(
    seq_result,
    chromosomes,
    samples,
    labels,
    mutation_filter,
    driver_mutations) {
  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_result) && ("mutations" %in% names(seq_result))) {

    # extract the field
    seq_res <- seq_result$mutations

    if (is.null(driver_mutations) && "parameters" %in% names(seq_result)) {
      if ("driver_mutations" %in% names(seq_result$parameters)) {
        driver_mutations <- seq_result$parameters$driver_mutations
      }
    }
  } else {
    seq_res <- seq_result
  }

  if ("sample" %in% names(seq_res)) {
    data <- seq_res
  } else {
    data <- seq_to_long(seq_res)
  }

  if (!is.null(labels)) {
    if (!is(labels, "data.frame")) {
      stop("The parameter \"labels\" must be a data frame when non-NULL.")
    }

    if (length(labels) != 1) {
      stop(paste0("The parameters \"labels\" must be a data frame ",
                  "with one column when non-NULL."))
    }

    if (nrow(labels) != nrow(seq_res)) {
      stop(paste0("The parameters \"seq_result\" and \"labels\"",
                  " must have the same number of rows."))
    }

    data["labels"] <- labels
    label_name <- names(labels)
  }

  if (!is.function(mutation_filter)) {
    stop("The parameter \"mutation_filter\" must be a function.")
  }
  data <- mutation_filter(data)

  chromosomes <- validate_chromosomes(seq_res, chromosomes)

  data <- data %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes)

  normal_data <- data %>% dplyr::filter(sample == "normal.sample")

  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample))) {
      wrong_names <- setdiff(samples, unique(data$sample))
      stop(paste("Invalid sample names in samples parameter:",
                 paste0(wrong_names, collapse = ", ")))
    }
    data <- data %>% dplyr::filter(.data$sample %in% samples)
  }

  result <- list("data" = data, "normal" = normal_data,
                 "driver_mutations" = driver_mutations)

  if (!is.null(labels)) {
    result["label_name"] <- label_name
  }

  return(result)
}

# setup data for VAF plotting
setup_VAF_data_for_plotting <- function(
    seq_result,
    chromosomes,
    samples,
    labels,
    mutation_filter,
    driver_mutations,
    cuts) {
  result <- setup_data_for_plotting(seq_result, chromosomes,
                                    samples, labels, mutation_filter,
                                    driver_mutations)

  data <- result$data %>%
    dplyr::filter(.data$VAF > 0) %>%
    dplyr::filter(.data$VAF <= max(cuts) & .data$VAF >= min(cuts))

  normal <- result$normal %>%
    dplyr::filter(.data$VAF > 0) %>%
    dplyr::filter(.data$VAF <= max(cuts) & .data$VAF >= min(cuts))

  filtered_result <- list("data" = data, "normal" = normal,
                          "driver_mutations" = result$driver_mutations)

  if (!is.null(labels)) {
    filtered_result["label_name"] <- result$label_name
  }

  return(filtered_result)
}

#' Plot a Variant Allele Frequency (VAF) histogram
#'
#' This function generates a histogram showing the distribution of Variant
#' Allele Frequency (VAF) across samples and chromosomes.
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
#' @param binwidth The width of the plot bins. When set to `NULL`, the
#'   function computes the most convenient bin width according to the
#'   maximum coverage reported in the data frame (default: `NULL`).
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
#' # plotting the VAF histogram without germinal mutations
#' plot_VAF_histogram(seq_results)
#'
#' # let us define a function to filter germinal and pre-neoplastic
#' # from the input data
#' library(dplyr)
#' filter_data <- function(data) {
#'   data %>% dplyr::filter(!classes %in% list("germinal",
#'                                             "pre-neoplastic"))
#' }
#'
#' # plotting the VAF histogram without germinal and pre-neoplastic
#' plot_VAF_histogram(seq_results, mutation_filter=filter_data)
#'
#' # plotting the VAF histogram filtering out VAFs below 0.02
#' plot_VAF_histogram(seq_results, cuts = c(0.02, 1))
#'
#' # plotting the VAF histogram with labels
#' plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
#'                    labels = seq_results$mutations["causes"])
#'
#' # avoid the driver mutation labels
#' plot_VAF_histogram(seq_results, cuts = c(0.02, 1),
#'                    labels = seq_results$mutations["causes"],
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
#' plot_VAF_histogram(f_seq, cuts = c(0.02, 1))
#'
#' # plotting the VAF histogram with labels
#' plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1))
#'
#' # use the driver codes in the driver mutation labels
#' plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1),
#'                    driver_mutations = phylo_forest$get_driver_mutations())
#'
#' # avoid the driver mutation labels
#' plot_VAF_histogram(f_seq, labels = f_seq["causes"], cuts = c(0.02, 1),
#'                    driver_mutation_labels = FALSE)
#'
#' # deleting the mutation engine directory
#' unlink('demo', recursive = T)
plot_VAF_histogram <- function(
  seq_result,
  chromosomes = NULL,
  samples = NULL,
  labels = NULL,
  binwidth = NULL,
  mutation_filter = filter_germinal,
  driver_mutations = NULL,
  driver_mutation_labels = TRUE,
  cuts = c(0, 1)
) {
  setup_data <- setup_VAF_data_for_plotting(seq_result, chromosomes, samples,
                                            labels, mutation_filter,
                                            driver_mutations, cuts)

  data <- setup_data$data

  if (!is.null(labels)) {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF,
                                             fill = labels)) +
      ggplot2::labs(col = NULL, fill = setup_data$label_name)
  } else {
    plot <- data %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF))
  }

  if (is.null(binwidth)) {
    binwidth <- freedman_diaconis_binwidth(data$VAF)
  }

  plot <- plot +
    ggplot2::geom_histogram(binwidth = binwidth, alpha = 0.5) +
    ggplot2::facet_grid(sample ~ chr, scales = "free_y") +
    ggplot2::coord_cartesian(ylim = c(0, NA), xlim = c(0, 1)) +
    ggplot2::ylab("count") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    ggplot2::theme(legend.position = "bottom")

  if (driver_mutation_labels) {
    plot <- add_driver_mutation_labels(plot, data, setup_data$driver_mutations,
                                       "VAF", 0)
  }

  return(plot)
}