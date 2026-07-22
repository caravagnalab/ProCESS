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

#' @name PhylogeneticForestNode$arising_mutations
#' @title Getting the mutations arising in the corresponding cell
#' @description This property is the data frame of the mutations arising
#'   the corresponding cell.
#' @return The data frame of the mutations arising the corresponding cell.
#' @example nobuild/roxygen/setups/PhylogeneticForestNode-setup.R
#' @examples
#'
#' # get the mutations arising in the corresponding cell
#' node$arising_mutations
#' @seealso <code>[PhylogeneticForestNode]</code>
NULL

#' @name PhylogeneticForestNode$get_genome
#' @title Getting the corresponding cell genome
#' @description This method computes the corresponding cell genome.
#' @details This method computes the corresponding cell genome by browsing
#'   the whole forest branch from the root down to the node. Whenever
#'   the genomes of many nodes are needed, using the node tour with genomes
#'   is preferable.
#' @return The genome of the corresponding cell.
#' @example nobuild/roxygen/setups/PhylogeneticForestNode-setup.R
#' @examples
#'
#' # the genome only has chromosome 22 because the forest was
#' # built by using the setup "demo"
#' node$get_genome()
#' @seealso <code>[PhylogeneticForestNode]</code>
NULL
