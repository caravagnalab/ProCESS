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

#' Get examples
#'
#' @description This function loads data structure examples.
#' @param name The name of the data structure example that is aimed.
#'   The supported examples are "SampleForest" and "PhylogeneticForest".
#' @return An example object of the specified data type.
#' @examples
#' # get an example of `SampleForest` object
#' forest <- ProCESS::example("SampleForest")
#'
#' # see the forest
#' forest
#'
#' # see the first nodes of the forest
#' head(forest$get_nodes())
#'
#' # load an example of `PhylogeneticForest` object
#' forest <- ProCESS::example("PhylogeneticForest")
#' @seealso `available_examples()`
#' @export
example <- function(name) {
  if (name == "SampleForest") {
    forest_path <- system.file("extdata", "sample_forest_example.sff",
                               package = "ProCESS")

    return(load_sample_forest(forest_path, quiet = TRUE))
  }
  if (name == "PhylogeneticForest") {
    forest_path <- system.file("extdata", "phylo_forest_example.pff",
                               package = "ProCESS")

    forest <- load_phylogenetic_forest(forest_path, quiet = TRUE)

    # we also need to set the reference path
    reference_path <- system.file("extdata", "example_ref.fasta",
                                  package = "ProCESS")

    forest$set_reference_path(reference_path)

    return(forest)
  }

  stop(paste0("Unknown example \"", name, "\"."))
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
  data.frame(
    "name" = c("SampleForest", "PhylogeneticForest"),
    "description" = c(paste0("The `SampleForest` object as build in ",
                             "[vignette::sampling",
                             "#two-populations-with-epigenetic-states]."),
                      paste0("The `PhylogeneticForest` object as build ",
                             "[vignette::mutations]."))
  )
}