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

#' @name SampleForest$get_samples_info
#' @title Retrieving the samples' information
#' @description This method retrieves information about
#'   the samples whose cells were used as leaves
#'   of the forest.
#' @return A data frame containing, for each sample collected
#'   during the simulation, the columns `name`, `time`, `id`,
#'   `ymin`, `xmin`, `ymax`, `xmax`, `tumour_cells`, and
#'   `tumour_cells_in_bbox`. The columns `ymin`, `xmin`, `ymax`,
#'   `xmax` report the boundaries of the sample bounding box, while
#'   `tumour_cells` and `tumour_cells_in_bbox` are the number of tumour
#'   cells in the sample and in the bounding box, respectively.
#' @examples
#' # use a forest example
#' forest <- example("SampleForest")
#'
#' # get information about the samples whose cells
#' # are the forest leaves
#' forest$get_samples_info()
#' @seealso [TissueSimulation$get_samples_info()],
#'   [PhylogeneticForest$get_samples_info()],
#'   <code>[SampleForest]</code>
NULL

#' @name SampleForest$get_species_info
#' @title Getting forest species
#' @description This method builds a data frame containing information
#'   about the simulated species.
#' @return A data frame reporting `mutant` and, if the simulation has
#'   epigenetic states, `epistate` for each registered species.
#' @examples
#' # use a forest example without epistates
#' forest <- example("SampleForest - no epistates")
#'
#' # get species information. Since the simulation has no epigenetic
#' # state, the species correspond to the mutants
#' forest$get_species_info()
#'
#' @examples
#' # use a forest example
#' forest <- example("SampleForest")
#'
#' # get species information
#' forest$get_species_info()
#' @seealso [MutationEngine$get_species_info()],
#'   [PhylogeneticForest$get_species_info()],
#'   <code>[SampleForest]</code>
NULL
