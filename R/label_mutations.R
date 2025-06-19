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

#' Label mutations using phylogenetic forest data
#'
#' @description
#' The function labels mutations using data about the cell in which it
#'   occurs for the first time.
#' @param seq_results The output of `simulate_seq()`.
#' @param phylo_forest The phylogenetic forest from which the sequencing
#'   was simulated.
#' @return A copy of the data frame `seq_results` added with the identifier of
#'   the cell in which the mutation occurs for the first time (column
#'   "`cell_id`"), the identifier of its ancestor (column "`ancestor`"), its
#'   mutant (column "`mutant`"), its epigenetic state (column "`epistate`"),
#'   its birth time (column "`birth_time`"), the sample that collected the
#'   cell whenever available (column "`sample`"), and a cell classification
#'   based on phylogenetic sticks (column "`label`").
#' @export
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' # simulate a tissue
#' sim <- SpatialSimulation()
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_size("A",1e4)
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B",1e4)
#'
#' # sample the tissue and build the sample forest
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#' forest = sim$get_samples_forest()
#'
#' # place the mutations
#' m_engine = MutationEngine(setup_code = "demo")
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 5e-8))
#' m_engine$add_mutant(mutant_name = "B",
#'                     passenger_rates = c(SNV = 5e-8))
#'
#' m_engine$add_exposure(time = 0, c(SBS1 = 0.2,SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 100, 10)
#'
#' # simulate sequencing without the normal sample
#' seq_results <- simulate_seq(phylo_forest, coverage = 100, write_SAM = F,
#'                             with_normal_sample = FALSE)
#'
#' library(dplyr)
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>% dplyr::filter(classes!="germinal")
#'
#' # label filtered mutations using phylogenetic forest data
#' labels <- label_mutations(f_seq, phylo_forest)
#' labels
#'
#' # plotting histogram of the VAF with phylogenetic labels
#' plot_VAF_histogram(f_seq, labels = labels["label"], cuts = c(0.02, 1))
label_mutations <- function(seq_results, phylo_forest) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_results) && ("mutations" %in% names(seq_results))) {

    # extract the field
    seq_res <- seq_results["mutations"]
  } else {
    seq_res <- seq_results
  }

  if (!inherits(phylo_forest, "Rcpp_PhylogeneticForest")) {
    stop("The second parameter must be a phylogenetic forest.")
  }

  cell_labels <- get_relevant_branches(phylo_forest)
  mutations_with_cell <- seq_res %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cell_id = phylo_forest$get_first_occurrences(Mutation(
      chr, chr_pos, ref, alt
    ))[[1]]) %>%
    dplyr::ungroup()
  labelled_mutations <- dplyr::full_join(mutations_with_cell,
                                         cell_labels,
                                         by = "cell_id") %>%
    dplyr::filter(!is.na(chr))
  return(labelled_mutations)
}