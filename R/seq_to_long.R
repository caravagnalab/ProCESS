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
#' "`ref`", "`alt`", "`NV`", "`DP`", "`VAF`", and "`sample_name`".
#' @export
#'
#' @examples
#' # Example data frame in wide format
#' seq_results <- data.frame(chr = c("chr1", "chr2"),
#'                           chr_pos = c(100, 200),
#'                           ref = c("A", "C"),
#'                           alt = c("T", "G"),
#'                           causes = c("SBS5", "SBS1"),
#'                           classes = c("germinal", "passenger"),
#'                           Sample.A.occurrences = c(10, 90),
#'                           Sample.A.coverage = c(100, 100),
#'                           Sample.A.VAF = c(0.1, 0.9),
#'                           normal.sample.occurrences = c(45, 52),
#'                           normal.sample.coverage = c(100, 100),
#'                           normal.sample.VAF = c(0.45, 0.52))
#' seq_results
#'
#' # Convert to long format
#' seq_to_long(seq_results)
seq_to_long <- function(seq_results) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_results) && ("mutations" %in% names(seq_results))) {

    # extract the field
    seq_res <- seq_results["mutations"]
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
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes",
            colnames(seq_res)[grepl(paste0(sn, "."),
                                    colnames(seq_res), fixed = TRUE)])

    # Rename columns and add sample_name column
    seq_res[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes",
                     "occurrences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)

  # Rename and reorder columns
  seq_df %>%
    dplyr::rename(from = chr_pos, DP = coverage,
                  NV = occurrences) %>%
    dplyr::mutate(to = from)
}