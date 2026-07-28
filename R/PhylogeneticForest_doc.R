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

#' @name PhylogeneticForest$partition_samples
#' @title Partitioning forest samples
#' @description This method partitions the samples in a phylogenetic forest.
#' @details This method partitions the samples in a phylogenetic forest
#'   according to a labelling function over the corresponding forest nodes.
#'   It works in-place by altering the phylogenetic forest from which the
#'   method is called.
#' @param labelling_function An R labelling function that maps any forest
#'   node to a labelling string.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # show information about samples
#' forest$get_samples_info()
#'
#' # the labelling function parameter has type `PhylogeneticForestNode`
#' birth_time_labelling <- function(node) {
#'   if (node$birth_time > 421) {
#'     return("YOUNG")
#'   }
#'
#'   if (node$birth_time > 321) {
#'     return("MIDDLE_AGED")
#'   }
#'
#'   return("OLD")
#' }
#'
#' # partition the samples according to the labelling function
#' forest$partition_samples(birth_time_labelling)
#'
#' # show the new samples
#' forest$get_samples_info()
#' @seealso <code>[PhylogeneticForestNode]</code>,
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_samples_info
#' @title Retrieving sample information
#' @description This method retrieves information about
#'   the samples whose cells were used as leaves
#'   of the sample forest.
#' @return A data frame containing, for each sample collected during the
#'   simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`,
#'    `ymax`, `ymax`, `xmax`, `tumour_cells`, `tumour_cells_in_bbox`,
#'   `DNA_quantity`, and `equivalent_normal_cells`. The columns `ymin`,
#'   `xmin`, `ymax`, and `xmax` report the boundaries of the sample
#'   bounding box, while `tumour_cells` and `tumour_cells_in_bbox` are the
#'   number of tumour cells in the sample and in the bounding box,
#'   respectively. Finally, `DNA_quantity` contains the overall number of
#'   tumoral bases, i.e., the sum of the lengths of all the alleles of all the
#'   sample tumoral cells, and `equivalent_normal_cells` contains the number
#'   of normal cells that contain as much DNA as the sample tumoral cells.
#'   The `DNA_quantity` is stored as a real number despite being a natural
#'   number as it usually exceeds the largest natural number that can be
#'   natively represented by R.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # show information about samples
#' forest$get_samples_info()
#' @seealso [SampleForest$get_samples_info()],
#'   [TissueSimulation$get_samples_info()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_driver_mutations
#' @title Getting the driver mutations
#' @description This method returns the applied driver mutations.
#' @return A data frame consisting in eight columns `order`, `type`,
#'   `CNA_type`, `chr`, `start`, `end`, `ref`, `alt`, `allele`, `allele_srd`,
#'   `cause`, and `code`. Each row in the data frame reports one driver
#'   mutations. The fields `cause` and `order` report the name of the mutant
#'   and the application order among the mutant driver mutations,
#'   respectively. The column `type` declares the mutation type and contains
#'   `SID`, `CNA`, or `WGD` when the mutation is an SNV/indel, a CNA, or
#'   a whole genome duplication, respectively. When the mutation is a CNA,
#'   the `CNA_type` can either be `A` (i.e., amplification) or `D`
#'   (i.e., deletion). When the mutation is not a WGD, the fields `chr`,
#'   `start`, and `end` contains the mutation chromosome, the initial and
#'   the final position on the chromosome, respectively. When the mutation
#'   is a SID , the fields `ref` and `alt` contains the mutation reference
#'   genome and alternate sequences, respectively. When the mutation is a
#'   SID or a CNA deletion, the field `allele` stores the allele in which
#'   the mutation was applied. When instead the mutation is a CNA
#'   amplification, the fields `allele` and `src_allele` reports the
#'   identifiers of the new allele and of the original allele, respectively.
#'   In all the remaining cases, the fields contains the value `NA`.
#'   Finally, the column `code` reports the mutation code.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # show information about samples
#' forest$get_driver_mutations()
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_species_info
#' @title Getting the species and their rates
#' @description This method returns the species and their rates.
#' @details This method returns the species and their rates during
#'   the simulation are returned in a data frame. The column `species`
#'   contains the species names; the columns `time`, `SNV_rate`,
#'   `indel_rate`, and `CNA_rate` store the time from which rates
#'   hold, and the corresponding the SNV, indel, and CNA rates,
#'   respectively.
#' @return A data frame reporting `species`, `time`, `SNV_rate`,
#'   `indel_rate`, and `CNA_rate` for each species.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # show information about samples
#' forest$get_species_info()
#' @seealso [MutationEngine$get_species_info()],
#'   [SampleForest$get_species_info()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_germline_subject
#' @title Getting the germline subject
#' @description This method returns a data frame reporting the germline
#'   subject name (column "sample"), population (column "pop"),
#'   super-population (column "super_pop"), and gender (column "gender").
#' @return The name of the subject whose germline is used.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the germline subject
#' forest$get_germline_subject()
#' @seealso [PhylogeneticForest$get_sampled_cell_mutations()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_sampled_cell_CNAs
#' @title Getting the sampled cells' CNAs
#' @description This method returns the CNAs of the sample cells.
#' @details This method builds a data frame representing all the CNAs
#'   in the cells sampled during the simulation and represented by
#'   the leaves of the phylogenetic forest.
#' @return A data frame reporting `cell_id`, `type` (`"A"` for amplifications
#'   and `"D"` for deletions), `chr`, `begin` (i.e., the first CNA
#'   locus in the chromosome), `end` (i.e., last CNA locus in the chromosome),
#'   `allele`, `src allele` (the allele origin for amplifications, `NA` for
#'   deletions), and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"`
#'   or `"pre-neoplastic"`).
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the sampled cell CNAs
#' CNAs <- forest$get_sampled_cell_CNAs()
#'
#' # print the first lines of the data frame
#' head(CNAs)
#' @seealso [PhylogeneticForest$get_sampled_cell_mutations()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_sampled_cell_mutations
#' @title Getting the sampled cells' mutations
#' @description This method returns the mutations of the sample cells.
#' @details This method builds a data frame representing all the SNV
#'   and the indel mutations in the cells sampled during the simulation
#'   and represented by the leaves of the phylogenetic forest.
#'   The data frame also reports the allele in which the mutations occur to
#'   support double occurrences due to CNAs.
#' @param with_germline A Boolean flag to report germline mutations too
#'   (default: `FALSE`).
#' @return A data frame reporting `cell_id`, `chr`, (i.e., the mutation
#'   chromosome), `from` (i.e., position in the chromosome), `allele`
#'   (in which the mutation occurs), `ref`, `alt`, `type` (i.e., either
#'   `"SNV"` or `"indel"`), `cause`, and `class` (i.e., `"driver"`,
#'   `"passenger"`, `"germinal"` or `"pre-neoplastic"`) for each mutation
#'   in the sampled cell genomes.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the sampled cell mutations
#' mutations <- forest$get_sampled_cell_mutations()
#'
#' # print the first lines of the data frame
#' head(mutations)
#' @seealso [PhylogeneticForest$get_sampled_cell_CNAs()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_germline_mutations
#' @title Getting the germinal mutations
#' @description This method returns the forest SNVs and indels.
#' @details Its builds a data frame representing all the germinal
#'   SNVs and indels of the cells represented in the phylogenetic forest.
#'   The data frame also reports the allele in which the mutations occur to
#'   support double occurrences due to CNAs.
#' @return A data frame reporting `chr`, `from` (i.e., the position in
#'   the chromosome), `allele` (in which the mutation occurs), `ref`,
#'   `alt`, `cause`, `type` (i.e., either `"SNV"` or `"indel"`) and
#'   `class` (i.e., `"germinal"`).
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the first germline mutations
#' head(forest$get_germline_mutations())
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_absolute_chromosome_positions
#' @title Getting the absolute chromosome positions
#' @description This method returns the absolute chromosome positions.
#' @details Its builds a data frame reporting the name, the length, and the
#'   initial and final absolute positions of each chromosome in the
#'   reference genome.
#' @return A data frame reporting the name (column `chr`), the length
#'   (column `length`), the initial absolute position (column `from`),
#'   and the final absolute position (column `to`) of each chromosome.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get absolute chromosome positions. Since this forest example was built
#' # by using one single chromosome, the resulting data frame contains only
#' # one line
#' forest$get_absolute_chromosome_positions()
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_exposures
#' @title Getting the timed exposure data frame
#' @description This method returns a data frame representing the exposure
#'   evolution over time.
#' @return A data frame reporting `time`, `signature`, `exposure` and,
#'   `type`.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the exposures used to build the forest
#' forest$get_exposures()
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_bulk_allelic_fragmentation
#' @title Getting the genome bulk allelic fragmentation
#' @description This method returns a data frame representing the bulk allelic
#'   fragmentation of the genome.
#' @param sample_name The name of the sample whose bulk allelic fragmentation
#'   is aimed (optional).
#' @return A data frame reporting, for each genomic fragment and for all
#'   the allelic type on the genomic fragment, the chromosome (`chr`),
#'   the first base position (`begin`), the last base position (`end`),
#'   the number of copy of the major and minor alleles (`major` and
#'   `minor`, respectively), and the ratio between the number of cells
#'   exhibiting this allelic type and the total number of cells in the
#'   sample.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the genome bulk allelic fragmentation
#' forest$get_bulk_allelic_fragmentation()
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_cell_allelic_fragmentation
#' @title Getting the cell allelic fragmentation data frame
#' @description This method returns a data frame representing the allelic
#'   fragmentation of each sampled cell.
#' @return A data frame reporting, for each cell, for each genomic fragment,
#'   and for all the allelic type on the genomic fragment, the cell
#'   identifier (`cell_id`), the chromosome (`chr`), the first base
#'   position (`begin`), the last base position (`end`), and the number
#'   of copy of the major and minor alleles (`major` and `minor`,
#'   respectively).
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # print the first rows of the cell allelic fragmentation
#' head(forest$get_cell_allelic_fragmentation())
#' @seealso `vignette("mutations")`,
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_first_occurrences
#' @title Getting the cell in which a mutation emerged
#' @description This method returns the identifier of the cell in which a
#'   mutation occurs for the first time.
#' @param mutation A mutation being a SNV, a indel, or a CNA.
#' @return The identifier of the cell in which a mutation
#'   occurs for the first time.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # consider a mutation in the forest
#' mutation <- Mutation("22", 35396109, ref = "CTGA", alt = "C")
#' mutation
#'
#' # get the identifiers of the cells in which the mutation arose, i.e.,
#' # they have the mutation, but their parent have not
#' cell_ids <- forest$get_first_occurrences(mutation)
#' cell_ids
#'
#' # get the corresponding node
#' node <- forest$get_node(cell_ids[[1]])
#'
#' # the mutation is among the arising mutation in the node
#' node$arising_mutations
#' @seealso `vignette("mutations")`,
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_mutation_statistics
#' @title Getting the statistics about mutations on each node
#' @description This method returns a data frame reporting the statistics about
#'   mutations on each node.
#' @return A data frame consisting of five columns `cell_id`, `new_SIDs`,
#'   `new_CNAs`, `total_SIDs`, and `total_CNAs`. Each row represents a
#'   node in the phylogenetic forest and reports the identifier of the
#'   corresponding cell and contains the number of mutations (`new_SIDs`) and
#'   CNAs (`new_CNAs`) appearing for the first time on the cell. Moreover, it
#'   show the total number of mutations and CNAs on the cell (`total_SIDs`
#'   and `total_CNAs`, respectively).
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the mutation statistics
#' forest$get_mutation_statistics()
#' @seealso <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$get_reference_path
#' @title Getting the reference genome path
#' @description This method returns the reference genome path.
#' @return The reference genome path.
#' @examples
#' # use a phylogenetic forest example
#' forest <- example("PhylogeneticForest")
#'
#' # get the reference path
#' forest$get_reference_path()
#' @seealso [PhylogeneticForest$set_reference_path()],
#'   <code>[PhylogeneticForest]</code>
NULL

#' @name PhylogeneticForest$set_reference_path
#' @title Setting the reference genome path
#' @description This method returns the reference genome path.
#' @return The reference genome path.
#' @seealso [PhylogeneticForest$get_reference_path()],
#'   <code>[PhylogeneticForest]</code>
NULL