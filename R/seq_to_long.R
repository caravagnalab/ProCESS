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

#' Convert sequencing results from wide to long format
#'
#' This function takes sequencing results in wide format and converts them into
#' a long format data frame. It extracts sample names from column names,
#' processes each sample separately, and then binds them together. Finally, it
#' renames and reorders columns to match the desired output format.
#'
#' @param seq_results A data frame containing sequencing results in wide format.
#' @return A data frame in long format with columns "`chr`", "`from`", "`to`",
#' "`ref`", "`alt`", "`NV`", "`DP`", "`VAF`", "`sample`", "`cause`", and
#' "`nature`".
#' @export
#'
#' @examples
#' # use a sequencing result example
#' seq_results <- example("Sequencing results")
#' head(seq_results)
#'
#' # Convert to long format
#' seq_long <- seq_to_long(seq_results)
#' head(seq_long)
seq_to_long <- function(seq_results) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_results) && ("mutations" %in% names(seq_results))) {

    # extract the field
    seq_res <- seq_results$mutations
  } else {
    seq_res <- seq_results
  }

  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_res)[grepl(".VAF",
                                                   colnames(seq_res),
                                                   fixed = TRUE)],
                           ".VAF") %>% unlist()

  # Process each sample separately to create a list of data frames
  seq_df <- lapply(sample_names, function(sn) {
    # Select relevant columns for the current sample
    cc <- c("chr", "from", "ref", "alt", "cause", "nature",
            colnames(seq_res)[grepl(paste0(sn, "."),
                                    colnames(seq_res), fixed = TRUE)])

    # Rename columns and add sample column
    seq_res[, cc] %>%
      `colnames<-`(c("chr", "from", "ref", "alt", "cause", "nature",
                     "NV", "DP", "VAF")) %>%
      dplyr::mutate(sample = sn)
  }) %>% do.call("bind_rows", .)

  # Rename and reorder columns
  seq_df %>%
    dplyr::rename(from = from, DP = DP,
                  NV = NV) %>%
    dplyr::mutate(to = from)
}