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

#' @name SampleForestNode$cell_id
#' @title Getting the identifier of the associated cell
#' @description This property stores the identifier of the
#'   cell associated to the node.
#' @return The identifier of the cell associated to the node.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the identifier of the cell associated to the node
#' node$cell_id
NULL

#' @name SampleForestNode$parent
#' @title Getting the parent node
#' @description This property stores the parent of the node.
#' @return The parent of the node.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the node's parent
#' node$parent
NULL

#' @name SampleForestNode$children
#' @title Getting the children of the node
#' @description This property stores the children of the node.
#' @return The children of the node.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the node's children
#' node$children
NULL

#' @name SampleForestNode$is_root
#' @title Checking whether the node is a root
#' @description This property is `TRUE` if and only if a root
#'    of the forest.
#' @return `TRUE` if and only if the node is a root of the forest.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # check whether the node is a root
#' node$is_root
#'
#' # get the node corresponding to the cell whose identifier is 1, i.e,
#' # the root
#' node <- forest$get_node(1)
#'
#' # check whether the node is a root
#' node$is_root
NULL

#' @name SampleForestNode$is_leaf
#' @title Checking whether the node is a leaf
#' @description This property is `TRUE` if and only if a leaf
#'    of the forest.
#' @return `TRUE` if and only if the node is a leaf of the forest.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # check whether the node is a leaf
#' node$is_leaf
#'
#' # get a tour over forest leaves
#' node_tour <- get_node_tour(forest, only_leaves = TRUE)
#'
#' # check whether the first node in the tour is a leaf
#' node_tour$node$is_leaf
NULL

#' @name SampleForestNode$sample_name
#' @title Getting the corresponding cell sample 
#' @description This property is the name of the sample that collected
#'   the corresponding cell.
#' @details This property is the name of the sample that collected
#'   the corresponding cell. If the node is not a leaf, it was not
#'   collected by any sample and this property is `NA`.
#' @return The name of of the sample that collected the corresponding
#'   cell if the corresponding cell was collected. Otherwise, `NA`.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # the corresponing cell was not collected and the property is `NA`
#' node$sample_name
#'
#' # get a tour over forest leaves
#' node_tour <- get_node_tour(forest, only_leaves = TRUE)
#'
#' # in this case, the node is a leaf and then was collected
#' node_tour$node$sample_name
NULL

#' @name SampleForestNode$birth_time
#' @title Getting the corresponding cell birth time
#' @description This property is the simulated time at which the corresponding
#'   time was born.
#' @return The simulated time at which the corresponding time was born.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the birth time
#' node$birth_time
NULL

#' @name SampleForestNode$death_time
#' @title Getting the corresponding cell death time
#' @description This property is the simulated time at which the corresponding
#'   time died.
#' @return The simulated time at which the corresponding time died or the 
#'   sampling time.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the death time
#' node$death_time
NULL

#' @name SampleForestNode$species_name
#' @title Getting the name of the corresponding cell species
#' @description This property is the name of the corresponding cell
#'   species.
#' @return The name of the corresponding cell species.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the name of the corresponding cell species
#' node$species_name
NULL

#' @name SampleForestNode$epistate_name
#' @title Getting the name of the corresponding cell epigenetic state
#' @description This property is the name of the corresponding cell
#'   species.
#' @return The name of the corresponding cell epigenetic state.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the name of the corresponding cell species
#' node$species_name
#'
#' # get the corresponding cell epigenetic state
#' node$epistate_name
NULL

#' @name SampleForestNode$mutant_name
#' @title Getting the name of the corresponding cell mutant
#' @description This property is the name of the corresponding cell
#'   mutant.
#' @return The name of the corresponding cell epigenetic state.
#' @example nobuild/roxygen/SampleForestNode-setup.R
#' @examples
#'
#' # get the name of the corresponding cell species
#' node$species_name
#'
#' # get the corresponding cell mutant
#' node$mutant_name
NULL
