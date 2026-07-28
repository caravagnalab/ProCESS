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

#' @name PhylogeneticForestNodeTour$genome
#' @title Getting the genome of the associated tour
#' @description This property stores the genome of the cell associated to
#'   the current node in the tour.
#' @details This property is optional and is available only if the
#'   `get_node_tour()`'s optional parameter `with_genomes` was set to `TRUE`.
#' @return The genome of the cell associated to the current node in the tour.
#' @example template_doc/setups/PhylogeneticForestNodeTour-setup.R
#' @examples
#'
#' #' # build a tour for the forest nodes
#' node_tour <- get_node_tour(forest, )
#'
#' # show the first node in the tour
#' node_tour$node
NULL
