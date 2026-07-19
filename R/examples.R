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

example_list <- list(
    "SampleForest" = list(
        type = "SampleForest",
        file = "p2_epi.sff",
        description = paste("A sample forest with two mutations,",
                            "two epigenetic states, and four samples.")),
    "PhylogeneticForest" = list(
        type = "PhylogeneticForest",
        file = "p2_epi.pff",
        description = paste("A phylogenetic forest with two mutations,",
                            "two epigenetic states, and four samples.")),
    "SampleForest - no epistates" = list(
        type = "SampleForest",
        file = "p2_no_epi.sff",
        description = paste("A sample forest with two mutations,",
                            "no epigenetic states, and four samples.")),
    "PhylogeneticForest - no epistates" = list(
        type = "PhylogeneticForest",
        file = "p2_no_epi.pff",
        description = paste("A phylogenetic forest with two mutations,",
                            "no epigenetic states, and four samples.")),
    "Sequencing results" = list(
        type = "sequencing results",
        file = "s_p2_epi.rds",
        description = paste("The result of a 10x sequencing of the",
                            "example \"PhylogeneticForest\"."))
)

#' Get examples
#'
#' @description This function loads data structure examples.
#' @param name The name of the data structure example that is aimed.
#'   The supported examples are "SampleForest" and "PhylogeneticForest".
#' @return An example object of the specified data type.
#' @examples
#' # get an example of `SampleForest` object
#' forest <- example("SampleForest")
#'
#' # see the forest
#' forest
#'
#' # see the first nodes of the forest
#' head(forest$get_nodes())
#'
#' # load an example of `PhylogeneticForest` object
#' forest <- example("PhylogeneticForest")
#' @seealso `available_examples()`
#' @export
example <- function(name) {
  index <- match(name, names(example_list))

  if (is.na(index)) {
    stop(paste0("Unknown example \"", name, "\". Check the available ",
                "examples using `available_examples()`."))
  }

  selected_example <- example_list[[index]]

  file_path <- system.file("extdata", selected_example$file,
                           package = "ProCESS")

  if (selected_example$type == "SampleForest") {
    return(load_sample_forest(file_path, quiet = TRUE))
  }

  if (selected_example$type == "PhylogeneticForest") {
    forest <- load_phylogenetic_forest(file_path, quiet = TRUE)

    # we also need to set the reference path
    reference_path <- system.file("extdata", "example_ref.fasta",
                                  package = "ProCESS")

    forest$set_reference_path(reference_path)

    return(forest)
  }

  if (selected_example$type == "sequencing results") {
    return(readRDS(file_path))
  }

  stop(paste0("Unknown example type \"", selected_example$type, "\"."))
}

#' Get the available data structure examples
#'
#' @description This function returns the available data structure
#'   examples.
#' @return A data frame describing the available examples and consisting of
#'   two columns: the names of the data structure examples (column `name`) and
#'   their descriptions (column `description`).
#' @examples
#' # get the data frame of the available data structure examples
#' available_examples()
#' @seealso `example()`
#' @export
available_examples <- function() {
  descriptions <- c()

  for (example in example_list) {
    descriptions <- c(descriptions, example$description)
  }

  data.frame(
    "name" = names(example_list),
    "description" = descriptions
  )
}