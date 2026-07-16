/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2026 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>

#include <Rcpp.h>

#include "mutation_engine.hpp"
#include "phylogenetic_forest.hpp"
#include "node_tour.hpp"
#include "sid.hpp"
#include "wg_doubling.hpp"

using namespace Rcpp;

RCPP_MODULE(Mutations)
{

//' @name Mutation
//' @title Either an SBS or an indel
    class_<SIDMut>("Mutation")
        .constructor()

//' @name Mutation$get_chromosome
//' @title Getting the mutation chromosome
//' @description This method identify the chromosome where the mutation
//'   occurred.
//' @return The identifier of the chromosome in which the mutation occurred.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the chromosome in which `snv` occurs (i.e., "X")
//' snv$get_chromosome()
        .method("get_chromosome", &SIDMut::get_chromosome,
                "Get the chromosome of the mutation")

//' @name Mutation$get_position_in_chromosome
//' @title Getting the mutation chromosome position
//' @description This method returns the position in the chromosome where
//'   the mutation occurred.
//' @return The position in chromosome where the mutation occurred.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the position in chromosome where `snv` occurs (i.e., 20002)
//' snv$get_position_in_chromosome()
        .method("get_position_in_chromosome", &SIDMut::get_position_in_chromosome,
                "Get the mutation position in the chromosome")

//' @name Mutation$get_ref
//' @title Getting the mutation reference sequence
//' @description This method returns the reference sequence that is
//'   altered by the mutation.
//' @return The reference sequence before the mutation.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the reference base in which `snv` occurs (i.e., "A")
//' snv$get_ref()
        .method("get_ref", &SIDMut::get_ref, "Get the mutation reference sequence")

//' @name Mutation$get_alt
//' @title Getting the mutation altered sequence
//' @description This method returns the sequence after the mutation occurred.
//' @return The sequence after the mutation occurred.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the sequence after `snv` occurs (i.e., "T")
//' snv$get_alt()
        .method("get_alt", &SIDMut::get_alt, "Get the mutation altered sequence")

//' @name Mutation$get_cause
//' @title Getting the mutation cause
//' @description This method returns the mutation cause.
//' @details Every mutation is associated to a cause depending on whether
//'   it is part of a genomic characterisation of a mutant or it is caused
//'   by a specific profile. This method returns such a cause whenever it is
//'   available.
//' @return The mutation cause.
//' @examples
//' # let us build a SNV without specifying any cause for it
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the cause of `snv` (i.e., "NA")
//' snv$get_cause()
//'
//' # we can also build a SNV, specifying a cause for it
//' snv <- SNV("X", 20002, "T", "A", cause = "SBS13")
//'
//' # get the cause of `snv` (i.e., "SBS13")
//' snv$get_cause()
        .method("get_cause", &SIDMut::get_cause, "Get the cause of the mutation")

//' @name Mutation$get_dataframe
//' @title Getting the mutation data frame
//' @description This method builds a data frame representing the mutation.
//' @details The data frame has the columns `chr`, `from`, `ref`, `alt`,
//'   `type` (i.e., `SNV` and `indel`), and `cause`.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' snv$get_dataframe()
        .method("get_dataframe", &SIDMut::get_dataframe,
                "Get a data frame representing the mutation")
        .method("show", &SIDMut::show);


//' @name SNV
//' @title Creating an SNV
//' @description This function creates SNVs.
//' @usage SNV(chr, chr_pos, alt, ref="?", allele = NULL, cause="")
//' @param chr The name of the chromosome in which the SNV occurs.
//' @param chr_pos The position in the chromosome where the SNV occurs.
//' @param alt The base after the mutation.
//' @param ref The base before the mutation (optional).
//' @param allele The allele in which the SNV must occur (optional).
//' @param cause The cause of the SNV (optional).
//' @seealso `Mutation()` for SNV and indel creation.
//' @examples
//' # create a SNV without specifying the cause and context
//' snv <- SNV("X", 20002, "T")
//' snv
//'
//' # create a SNV and do not specify the cause
//' snv <- SNV("X", 20002, "T", "A")
//' snv
//'
//' # create a SNV that must be place in allele 1
//' snv <- SNV("X", 20002, "T", allele = 1)
//' snv
//'
//' # create a SNV with a cause
//' snv <- SNV("X", 20002, "T", cause = "SBS1")
//' snv
    function("SNV", &SIDMut::build_SNV,
             List::create(_["chr"], _["from"], _["alt"], _["ref"] = "?",
                          _["allele"] = R_NilValue, _["cause"] = ""),
             "Create a single nucleotide variation (SNV)");

//' @name Mutation
//' @title Creating a SNV or a indel
//' @description This function creates SNVs and indels.
//' @details This function generalizes the function `SNV()` by constructing
//'   SNVs and indels. However, it necessitates the specification of the
//'   reference sequence, whereas `SNV()` can infer it from the reference
//'   sequence itself.
//'
//'   Another distinction between this function and `SNV()` lies in the order
//'   of the `ref`-`alt` parameter: in `SNV()`, the alt parameter precedes the
//'   optional ref parameter, while `Mutation()` adopts the reverse order.
//' @param chr The name of the chromosome in which the indel occurs.
//' @param from The position in the chromosome where the indel occurs.
//' @param ref The reference sequence.
//' @param alt The mutation altered sequence.
//' @param allele The allele in which the mutation must occur (optional).
//' @param cause The cause of the mutation (optional).
//' @seealso `SNV()` for SNV creation.
//' @examples
//' # create a deletion without specifying the cause
//' mutation <- Mutation("X", 20002, "TAC", "T")
//' mutation
//'
//' # create an insertion and do not specify the cause
//' mutation <- Mutation("X", 20002, "A", "AT")
//' mutation
//'
//' # create an insertion that must be place in allele 1
//' mutation <- Mutation("X", 20002, "A", "AT", allele = 1)
//' mutation
//'
//' # create an insertion with a cause
//' mutation <- Mutation("X", 20002, "A", "AT", cause = "SBS1")
//' mutation
    function("Mutation", &SIDMut::build_SID,
             List::create(_["chr"], _["from"], _["ref"], _["alt"],
                          _["allele"] = R_NilValue, _["cause"] = ""),
             "Create either an SNV or a indel");

//' @name CNA
//' @title Creating a CNA
//' @description This function creates a CNA.
//' @usage CNA(type, chr, chr_pos, len, allele = NULL, src_allele = NULL)
//' @param type The CNA type: either `"A"` or `"D"` for amplification and
//'   deletion, respectively.
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param from The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the CNA occurs. (optional)
//' @param src_allele The allele from which the region is amplified. (optional,
//'   for amplification only)
//' @seealso  [Amplification()] to build an amplification;
//'   [Deletion()] to build a deletion.
//' @examples
//' # create an amplification
//' cna <- CNA("A", "X", 20002, 100)
//'
//' cna
//'
//' # create a deletion from the allele 0
//' cna <- CNA("D", "Y", 101310, 205, allele = 0)
//'
//' cna
    function("CNA", &CNA::build_CNA,
             List::create(_["type"], _["chr"], _["from"], _["len"],
                          _["allele"] = R_NilValue, _["src_allele"] = R_NilValue),
             "Create a copy number alteration (CNA)");

//' @name Amplification
//' @title Creating a CNA amplification
//' @description This function creates a CNA amplification.
//' @usage Amplification(chr, chr_pos, len, allele = NULL, src_allele = NULL)
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param from The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the amplification is placed. (optional)
//' @param src_allele The allele from which the region is amplified. (optional)
//' @seealso  [Deletion()] to build a deletion; [CNA()]
//'   to build both amplifications and deletions.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100)
//'
//' cna
    function("Amplification", &CNA::build_amplification,
             List::create(_["chr"], _["from"], _["len"], _["allele"] = R_NilValue,
                          _["src_allele"] = R_NilValue),
             "Create a CNA amplification");

//' @name Deletion
//' @title Creating a CNA deletion
//' @description This function creates a CNA deletion.
//' @usage Deletion(chr, chr_pos, len, allele = NULL)
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param from The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the deletion occurs. (optional)
//' @seealso  [Amplification()] to build an amplification;
//'   [CNA()] to build both amplifications and deletions.
//' @examples
//' # create a deletion CNA
//' cna <- Deletion("Y", 40020, 200)
//'
//' cna
    function("Deletion", &CNA::build_deletion,
             List::create(_["chr"], _["from"], _["len"], _["allele"] = R_NilValue),
             "Create a CNA deletion");

//' @name CNA
//' @title A copy number alteration
    class_<CNA>("CNA")
        .constructor()

//' @name CNA$get_chromosome
//' @title Getting the CNA chromosome
//' @description This method returns the identifier of the chromosome
//'   where the CNA occurred.
//' @return The identifier of the chromosome in which the CNA occurred.
//' @examples
//' # create an amplification CNA
//' cna <- CNA("A", "X", 20002, 100)
//'
//' # get the chromosome in which `cna` occurs (i.e., "X")
//' cna$get_chromosome()
        .method("get_chromosome", &CNA::get_chromosome, "Get the chromosome of the CNA")

//' @name CNA$get_position_in_chromosome
//' @title Getting the CNA chromosome position
//' @description This method returns the position in chromosome
//'   where the CNA occurred.
//' @return The position in chromosome where the CNA occurred.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the position in chromosome where `cna` occurs (i.e., 20002)
//' cna$get_position_in_chromosome()
        .method("get_position_in_chromosome", &CNA::get_position_in_chromosome,
                "Get the CNA position in the chromosome")

//' @name CNA$get_length
//' @title Getting the CNA length
//' @description This method returns the CNA length.
//' @return The CNA length.
//' @examples
//' # create an amplification CNA
//' cna <- CNA("A", "X", 20002, 100)
//'
//' # get the length of `cna` (i.e., 100)
//' cna$get_length()
        .method("get_length", &CNA::get_length, "Get the CNA length")

//' @name CNA$get_allele
//' @title Getting the CNA allele
//' @description This method returns the identifier of the allele in
//'    which the CNA occurred.
//' @details If the CNA is an amplification corresponds to the new
//'    allele identifier. If, instead, the CNA is a deletion is the
//'    identifier of the allele on which the deletion occurred.
//' @return The allele in which CNA occurred.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the allele in which `cna` occurs (i.e., 1)
//' cna$get_allele()
        .method("get_allele", &CNA::get_allele, "Get the alteration allele")

//' @name CNA$get_src_allele
//' @title Getting the CNA source allele
//' @description This method returns the identifier of the allele from
//'    which the CNA is copied.
//' @return The allele from which CNA is copied.
//' @examples
//' # create an amplification CNA
//' amp_cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get allele from which `amp_cna` is copied (i.e., 0)
//' amp_cna$get_src_allele()
//'
//' # create a deletion CNA
//' del_cna <- Deletion("Y", 40020, 200, 0)
//'
//' # the deletions have no sources and the method returns NA
//' del_cna$get_src_allele()
        .method("get_src_allele", &CNA::get_src_allele,
                "Get the source allele (for amplifications)")

//' @name CNA$get_dataframe
//' @title Getting the CNA data frame
//' @description This method builds a data frame representing the CNA.
//' @details The data frame contains the  columns `chr`, `from`,
//'   `length`, `alt_base`, `allele`", `src.allele`, and `type`.
//' @examples
//' # create an amplification CNA
//' amp_cna <- Amplification("X", 20002, 100)
//'
//' amp_cna$get_dataframe()
//'
//' # create a deletion CNA
//' del_cna <- Deletion("Y", 40020, 200, 0)
//'
//' del_cna$get_dataframe()
        .method("get_dataframe", &CNA::get_dataframe,
                "Get a data frame representing the CNA")
        .method("show", &CNA::show);

//' @name WholeGenomeDoubling
//' @title Whole genome doubling events
//' @description A whole genome doubling event (WGD)
//'   produces the simultaneous duplication of all the
//'   chromosome allele in a genome.
//' @keywords internal
    class_<WholeGenomeDoubling>("WholeGenomeDoubling")
        .constructor()
        .method("show", &WholeGenomeDoubling::show);

//' @name MutationEngine
//' @title Generating phylogenetic forests
//' @description A mutation engine can label every node of a samples
//'   forest by mutations and produce a consistent phylogenetic forest.
//' @details The mutations are randomly generated according to three
//'   factors:
//'   1. the mutational rates of the species involved in the samples
//'   forest
//'   2. the genotypical specification of the mutants involved in the
//'   sample forest, i.e., the somatic mutations characterising
//'   the mutant genotypes
//'   3. the SBS and ID signatures active along the species simulation
//'
//'   The data of points 1 and 2 are provided to the mutation engine by
//'   the method [MutationEngine$add_mutant()].
//'   Instead, the active signatures are defined by using the method
//'   [MutationEngine$add_exposure()].
//'
//'   The initialisation of a `MutationEngine` object requires a reference
//'   sequence and the SBS and ID mutational signatures. An SBS index and
//'   a ID index of the reference sequence are then automatically built.
//'   This process may take time depending on the size of the reference
//'   sequence. Hence, the downloaded files together with the built indices
//'   are saved on the disk for subsequent `MutationEngine` constructions.
//'
    class_<MutationEngine>("MutationEngine")

//' @name MutationEngine$infinite_sites_model
//' @title Switching on and off the infinite sites model
//' @description This property enables/disables the infinite sites model.
//' @details When it is `TRUE`, the infinite sites model is enabled and
//'   new mutations are exclusively placed in mutation-free loci.
//' @examples
//' # create a demonstrative mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # the infinite sites model is enabled by default
//' m_engine$infinite_sites_model
//'
//' # the infinite sites model can be disabled
//' m_engine$infinite_sites_model <- FALSE
//'
//' m_engine$infinite_sites_model
        .property("infinite_sites_model", &MutationEngine::get_infinite_sites_model,
                  &MutationEngine::set_infinite_sites_model,
                  "A flag to enable/disable the infinite sites model")

//' @name MutationEngine$add_exposure
//' @title Adding an exposure to the mutation engine
//' @description This method adds an exposure to the mutation engine.
//' @details The exposure will be used to establish the probability
//'   for a passenger mutation to occur depending on its context.
//'
//'   Each exposure is associated to a time that is the simulated
//'   time in which the set is adopted.
//'   If a time is provided the exposure is used from the specified
//'   time on up to the successive exposure change. When an exposure
//'   is added to the mutation engine without specifying the time,
//'   its time is 0.
//' @param time The simulated time at which the exposure is adopted
//'   (optional).
//' @param exposure An exposure for the specified mutation type, i.e.,
//'   a discrete probability distribution over a set of signature.
//'   The indel and SNV exposures can be specified in the same list.
//' @examples
//' # create a demonstrative mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add a default set of coefficients that will be used from simulated
//' # time 0 up to the successive coefficient change. The indel and SNV
//' # exposures can be specified in the same list.
//' m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7, ID2 = 0.2, ID3 = 0.3,
//'   ID20 = 0.5))
//'
//' # add a default set of coefficients that will be used from simulated
//' # time 3.2 up to the end of the simulation.
//' m_engine$add_exposure(3.2, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))
//'
//' m_engine
        .method("add_exposure",
                (void (MutationEngine::*)(const List &))(&MutationEngine::add_exposure),
                "Add an exposure")
        .method("add_exposure",
                (void (MutationEngine::*)(const double &, const List &))(
                    &MutationEngine::add_exposure),
                "Add an exposure")

//' @name MutationEngine$add_mutant
//' @title Adding a mutant specification
//' @description This method adds a mutant specification to the mutation engine.
//' @details The users must use it to specify the name and the genomic
//'   characterisation (i.e., SNVs, indels, CNAs, and whole genome doubling
//'   events (WGD)) of all the simulated mutants together with the mutation
//'   rates of its species.
//'   The driver mutations are applied to the mutant progenitor's genome
//'   respecting the specification order.
//' @param mutant_name The mutant name.
//' @param passenger_rates The list of the passenger rates whose names are the
//'   epigenetic states of the species or a single rate, if the mutant
//'   does not have an epigenetic state.
//' @param drivers The list of the driver SNVs, indels, CNAs, and the whole
//'   genome doubling events (WGD) characterizing the mutant (optional).
//' @seealso [MutationEngine$change_rates_from()]
//' @examples
//' # create a demonstrative mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # define a list of mutations
//' d_mutations <- list("DGCR8 P26L",
//'                     Mutation("22", 16085675, "GCCTCCCGA", "G"),
//'                     "EP300 S2346del",
//'                     WGD,
//'                     CNA(type = "A", chr = "22", from = 10303470,
//'                         len = 200000),
//'                     SNV("22", 23657587, "C"),
//'                     CNA("D", "22", 5010000, 200000))
//'
//' # add the mutant "A" characterized by the mutations in `d_mutations`. The
//' # mutations are applied according to `d_mutations`'s order. The mutant has
//' # one epigenetic states and its species "A[E1]" and "A[E2]" have passenger
//' # SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
//' # 1e-11, respectively.
//' m_engine$add_mutant("A", passenger_rates =c(SNV = 1e-9, indel = 1e-10),
//'                     drivers = d_mutations)
//'
//' m_engine
//'
//' # add the mutant "B" characterized by the mutations in `d_mutations`. The
//' # mutations are applied according to `d_mutations`'s order. The mutant has
//' # two epigenetic states and its species "B[E1]" and "B[E2]" have passenger
//' # SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
//' # 1e-11, respectively.
//' m_engine$add_mutant("B", list("E1" = c(SNV = 1e-9, indel = 1e-10),
//'                               "E2" = c(SNV = 3e-8, CNA = 1e-11)),
//'                     drivers = d_mutations)
//'
//' m_engine
        .method("add_mutant",
                (void (MutationEngine::*)(const std::string &,
                                          const Rcpp::List &))(
                    &MutationEngine::add_mutant),
                "Add mutant")
        .method("add_mutant",
                (void (MutationEngine::*)(const std::string &,
                                          const Rcpp::List &,const Rcpp::List &))(
                    &MutationEngine::add_mutant),
                "Add mutant")

//' @name MutationEngine$change_rates_from
//' @title Change the passenger rates from a specified time
//' @description This method changes the passenger rates from a specified time.
//' @details This method changes the passenger rates from a specified time. The
//'   rates before the specified time and those of the unspecified epigenetic
//'   states are not affected.
//' @param time The time from which the passenger rates are set.
//' @param mutant_name The mutant name.
//' @param passenger_rates The list of the passenger rates whose names are the
//'   epigenetic states of the species or a single rate, if the mutant
//'   does not have an epigenetic state.
//' @seealso [MutationEngine$add_mutant()]
//' @examples
//' # create a demonstrative mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # define a list of mutations
//' d_mutations <- list("DGCR8 P26L",
//'                     Mutation("22", 16085675, "GCCTCCCGA", "G"),
//'                     "EP300 S2346del",
//'                     WGD,
//'                     CNA(type = "A", chr = "22", from = 10303470,
//'                         len = 200000),
//'                     SNV("22", 23657587, "C"),
//'                     CNA("D", "22", 5010000, 200000))
//'
//' # add the mutant "A" characterized by the mutations in `d_mutations`. The
//' # mutations are applied according to `d_mutations`'s order. The mutant has
//' # two epigenetic states and its species "A[E1]" and "A[E2]" have passenger
//' # SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
//' # 1e-11, respectively.
//' m_engine$add_mutant("A", list("E1" = c(SNV = 1e-9, indel = 1e-10),
//'                               "E2" = c(SNV = 3e-8, CNA = 1e-11)),
//'                     drivers = d_mutations)
//'
//' m_engine
//'
//' # change the rates of "A[E1]" from time 10
//' m_engine$change_rates_from(10, "A", list("E1" = c(SNV = 2e-9, indel = 4e-9)))
//'
//' # ... and those of "A[E2]" from time 13
//' m_engine$change_rates_from(13, "A", list("E2" = c(SNV = 2e-9)))
//'
//' m_engine
        .method("change_rates_from",
                (void (MutationEngine::*)(const CLONES::Time, const std::string &,
                                          const Rcpp::List &))(
                    &MutationEngine::change_rates_from),
                "Change passenger rates from a given timestamp")

//' @name MutationEngine$place_mutations
//' @title Placing the mutations
//' @description This methods places mutations on a sample forest.
//' @details Each node of a sample forest is labelled by the
//'   mutations occurring in the cell represented by the node itself
//'   and produces a phylogenetic forest.
//' @param sample_forest A sample forest.
//' @param num_of_pre_neoplastic_SNVs The number of pre-neoplastic SNVs.
//' @param pre_neoplastic_SNV_signature_name The name of the SNV signature
//'   for the pre-neoplastic SNV generation (optional).
//' @param num_of_pre_neoplastic_indels The number of pre-neoplastic indels.
//' @param pre_neoplastic_indel_signature_name The name of the indel signature
//'   for the pre_neoplastic indel generation.
//' @param seed The seed for random number generator (optional).
//' @return A phylogenetic forest whose structure corresponds to
//'   `sample_forest`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # build a tissue simulation
//' sim <- TissueSimulation()
//'
//' # add mutant "A" and set its rates
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.01))
//'
//' # place a cell of species "A" in position (500,500)
//' sim$place_cell("A", 500, 500)
//'
//' # run the simulation until "A" accounts for less than 50000 cells
//' sim$run_up_to_size("A", 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner = c(450, 475),
//'                  upper_corner = c(500, 550))
//'
//' # build the sample forest
//' sample_forest <- sim$get_sample_forest()
//'
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add the mutant "A" to the engine
//' m_engine$add_mutant("A", c(SNV = 3e-9), list(SNV("22", 12028576, "G")))
//'
//' # add the default set of SNV signature coefficients
//' m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7, ID2 = 0.3, ID21 = 0.5,
//'                         ID3 = 0.2))
//'
//' # place the mutations on the sample forest assuming 1000 pre-neoplastic
//' # SNVs and 500 indels
//' phylogenetic_forest <- m_engine$place_mutations(sample_forest, 1000, 500)
//'
//' phylogenetic_forest
        .method("place_mutations",
                (PhylogeneticForest (MutationEngine::*)(
                    const SampleForest &forest, const size_t &num_of_pre_neoplastic_SNVs,
                    const size_t &num_of_pre_neoplastic_indels))(
                    &MutationEngine::place_mutations),
                "Placing mutations on a SampleForest")
        .method("place_mutations",
                (PhylogeneticForest (MutationEngine::*)(
                    const SampleForest &forest, const size_t &num_of_pre_neoplastic_SNVs,
                    const size_t &num_of_pre_neoplastic_indels,
                    const SEXP &seed))(&MutationEngine::place_mutations),
                "Placing mutations on a SampleForest")
        .method("place_mutations",
                (PhylogeneticForest (MutationEngine::*)(
                    const SampleForest &forest, const size_t &num_of_pre_neoplastic_SNVs,
                    const std::string &pre_neoplastic_SNV_signature_name,
                    const size_t &num_of_pre_neoplastic_indels,
                    const std::string &pre_neoplastic_indel_signature_name))(
                    &MutationEngine::place_mutations),
                "Placing mutations on a SampleForest")
        .method("place_mutations",
                (PhylogeneticForest (MutationEngine::*)(
                    const SampleForest &forest, const size_t &num_of_pre_neoplastic_SNVs,
                    const std::string &pre_neoplastic_SNV_signature_name,
                    const size_t &num_of_pre_neoplastic_indels,
                    const std::string &pre_neoplastic_indel_signature_name,
                    const SEXP &seed))(&MutationEngine::place_mutations),
                "Placing mutations on a SampleForest")

//' @name MutationEngine$get_genome_info
//' @title Getting the genome information
//' @description This method returns information about the genome.
//' @details This method returns a data frame reporting the name
//'    (column `name`), the size (column `size`), and the number
//'    of alleles (column `num_of_alleles`) of each chromosome.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the genome information
//' m_engine$get_genome_info()
        .method("get_genome_info", &MutationEngine::get_genome_info)

//' @name MutationEngine$get_active_germline
//' @title Getting the active germline subject
//' @description This method returns the active germline subject.
//' @details The active germline subject is returned as a
//'   data frame in which the column `sample` reports the
//'   subject name, the columns `pop` and `super_pop` contain the
//'   subject population and super population, respectively, and
//'   the column `gender` declares the subject gender.
//' @return A data frame the active germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to
//'   get the available germline subjects;
//'   [MutationEngine$set_germline_subject()] to set the
//'   active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_active_germline(), 5)
        .method("get_active_germline", &MutationEngine::get_active_germline)

//' @name MutationEngine$set_germline_subject
//' @title Setting the germline subject
//' @description This method sets the germline subject.
//' @details The subject must be one among those reported by
//'   [MutationEngine$get_germline_subjects()].
//' @return Set the germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to
//'   get the available germline subjects;
//'   [MutationEngine$get_active_germline()] to get the
//'   active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # set the active germline subject data frame
//' m_engine$set_germline_subject("NA18941")
        .method("set_germline_subject", &MutationEngine::set_germline_subject)

//' @name MutationEngine$get_germline_subjects
//' @title Getting the germline subjects
//' @description This method returns the available germline subjects.
//' @details The germline subjects method returns a data frame
//'   containing the available germline subjects. The column `sample`
//'   reports the subject name; the columns `pop` and `super_pop`
//'   contain the subject population and super population,
//'   respectively; the column `gender` declares the subject gender.
//' @return A data frame the available germline subjects.
//' @seealso [MutationEngine$get_active_germline()] to get the
//'   available germline subjects;
//'   [MutationEngine$set_germline_subject()]
//'   to set the active germline.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_germline_subjects(), 5)
        .method("get_germline_subjects", &MutationEngine::get_germline_subjects)

//' @name MutationEngine$get_population_descriptions
//' @title Getting the population descriptions
//' @description This method returns the population descriptions.
//' @details The population descriptions are stored in a
//'   data frame describing the populations. The column `code`
//'   contains the population codes; the columns `description`
//'   and `long description` report a brief and a long
//'   description for the populations, respectively.
//' @return A data frame containing the population descriptions.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_population_descriptions(), 5)
        .method("get_population_descriptions",
                &MutationEngine::get_population_descriptions)

//' @name MutationEngine$get_species_info
//' @title Getting the registered species and their rates
//' @description This method returns the registered species and
//'   their rates.
//' @details The registered species and their rates during the
//'   simulation are returned in a data frame. The column
//'   `species` contains the species names; the columns `time`,
//'   `SNV_rate`, `indel_rate`, and `CNA_rate` store the time
//'   from which rates hold, and the corresponding the SNV,
//'   indel, and CNA rates, respectively.
//' @return A data frame containing the registered species rates.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_species_info(), 5)
//' @seealso [PhylogeneticForest$get_species_info()]
        .method("get_species_info", &MutationEngine::get_species_info)

//' @name MutationEngine$get_SNV_signatures
//' @title Getting the SNV signatures
//' @description This method returns the available SNV
//'   signatures.
//' @details The signatures are returned in a data frame
//'   containing the available SNV signatures and the
//'   corresponding mutation probability. The first column
//'   ("Type") describes a mutation in a context, while each
//'   of the remaining columns contains the probabilities
//'   of the mutations for one of the available SNV
//'   signatures.
//' @return A data frame containing the available SNV
//'   signatures.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the indel data frame
//' head(m_engine$get_SNV_signatures(), 5)
        .method("get_SNV_signatures", &MutationEngine::get_SNV_signatures_dataframe,
                "Get the SNV signatures data frame")

//' @name MutationEngine$get_indel_signatures
//' @title Getting the indel signatures
//' @description This method returns the available indel
//'   signatures.
//' @details The signatures are returned in a data frame
//'   containing the available indel signatures together with
//'   the corresponding mutation probability. The first column
//'   ("Type") describes a mutation in a context, while each
//'   of the remaining columns contains the probabilities of
//'   the mutations for one of the available indel signatures.
//' @return A data frame containing the available  indel
//'   signatures.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the indel data frame
//' head(m_engine$get_indel_signatures(), 5)
        .method("get_indel_signatures", &MutationEngine::get_indel_signatures_dataframe,
                "Get the indel signatures data frame")

//' @name MutationEngine$get_known_drivers
//' @title Getting the known driver mutations
//' @description This method returns the known driver
//'   mutations.
//' @details The mutation are returned in a data frame reporting
//'   the known driver mutations together with their types,
//'   associated tumours, affected genes, and code name. The
//'   first three columns ("`chr`", `from`, and `to`)
//'   report the mutation chromosome, the initial position
//'   and the final position, respectively. The next three
//'   columns ("`ref`", `alt`, and `mutation_type`)
//'   describe the reference sequence, the altered sequence,
//'   and the type of the mutation. The last four columns
//'   ("`driver_gene`", `driver_code`, `driver_CDS`, and
//'   `tumour_type`) detail the affected gene, the driver
//'   code, which can be used to specify the mutation when
//'   adding a mutant to the mutation engine, the variant code,
//'   and the tumour type associated to the mutation.
//' @return A data frame containing the known driver.
//' @seealso [MutationEngine$add_mutant()]
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the known driver data frame
//' head(m_engine$get_known_drivers(), 5)
        .method("get_known_drivers", &MutationEngine::get_known_driver_mutations,
                "Get the known driver data frame")

        .method("show", &MutationEngine::show);

//' @name MutationEngine
//' @title Creating a mutation engine
//' @description This function downloads and sets up the data
//'   requires by a mutation engine. Finally, it builds mutation
//'   engine itself.
//'
//' @details There are two building modalities: the first one is more
//'   general, but it requires to specify all the data sources; the
//'   second one adopts some pre-set configurations, but it is
//'   sufficient in many cases.
//'
//'   The first building modality requires to specify the directory in
//'   which the data must be saved, the path or URL of the reference
//'   sequence, the mutational signatures, the driver SNVs file, the
//'   passenger CNAs file, and the germline data directory through the
//'   `directory`, `reference_src`, `SBS_src`, `drivers_src`,
//'   `passenger_CNAs_src`, and `germline_src`, respectively.
//'
//'   The second building modality exclusively requires a set-up code
//'   (parameter `setup_code`). The list of supported set-up codes can
//'   be obtained by using the function
//'   [get_mutation_engine_codes()].
//'
//'   Whenever the mutational signatures are meant to be downloaded from
//'   the COSMIC site, a valid COSMIC account is needed and can be
//'   provided by the parameter `COSMIC_account`.
//'
//'   The number of context sampling is an optional parameter that allows
//'   sampling the reference contexts while building the context index.
//'   This parameter, which is set to 100 by default, specifies how many
//'   occurrences of the same context must be identified before adding
//'   one of them to the context index. The larger the number of context
//'   sampling, the larger the context index. On the other side, the
//'   lower the number of context sampling, the lower the number of sites
//'   in the reference genome that can be affected by simulated
//'   mutations.
//'
//'   If the parameters of a mutation engine construction match those
//'   of a previous construction, then the corresponding reference
//'   sequence, the SBS file, and the previously built context index
//'   are loaded from the set-up directory avoiding further
//'   computations.
//' @usage
//' MutationEngine(directory, reference_src,
//'                SBS_signatures_src, indel_signatures_src,
//'                drivers_src, passenger_CNAs_src,
//'                germline_src)
//'
//' MutationEngine(setup_code="demo")
//' @export
//' @param setup_code The set-up code (optional).
//' @param directory The set-up directory (mandatory when `setup_code` is
//'   *not* provided).
//' @param reference_src The reference genome path or URL (mandatory when
//'   `setup_code` is *not* provided).
//' @param SBS_signatures_src The SBS signature file path or URL (mandatory
//'   when `setup_code` is *not* provided).
//' @param indel_signatures_src The indel signature file path or URL (mandatory
//'   when `setup_code` is *not* provided).
//' @param drivers_src The driver mutation file path or URL (mandatory
//'   when `setup_code` is *not* provided).
//' @param passenger_CNAs_src The passenger CNAs file path or URL (mandatory
//'   when `setup_code` is *not* provided).
//' @param germline_src The germline directory path or URL (mandatory when
//'   `setup_code` is *not* provided).
//' @param germline_subject The germline subject (optional).
//' @param context_sampling The number of reference contexts per context in
//'   the index (optional: default value is 30).
//' @param COSMIC_account A named list containing "email" and "password" of
//'   a valid COSMIC account (required to download mutational signatures
//'   from COSMIC site).
//' @param max_index_size The maximum size of an admitted indel and, as a
//'   consequence, the maximum size of a motif stored in the repeated
//'   sequence index (optional: default value is 50).
//' @param max_repetition_storage The maximum number of repetitions per type
//'   stored in the repeated sequence index (optional: default value is
//'   500000).
//' @param driver_CNA_min_distance The minimum distance between a driver
//'   mutation and a passenger CNA.
//' @param tumour_type The type of tumour. This is currently used to select
//'   the admissible passenger CNAs. If any passenger CNA in the dataset is
//'   admissible, use the the empty string `""` (optional: default value is
//'   `""`).
//' @param avoid_homozygous_losses An optional Boolean flag to avoid
//'   homozygous losses. When set to `TRUE`, passenger CNAs will be
//'   exclusively applied to regions covered by two alleles at least.
//'   (default: TRUE).
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso [get_mutation_engine_codes()] provides a list of
//'   the supported set-up codes;
//'   [MutationEngine$get_germline_subjects()] to get the
//'   available germline subjects;
//'   [MutationEngine$set_germline_subject()] to set the active
//'   germline subject; [MutationEngine$get_active_germline()] to
//'   get the active germline subject.
//' @return A mutation engine object.
//' @examples
//' # set the reference and SBS URLs
//' reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/release-111/",
//'   "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
//'   "dna.chromosome.22.fa.gz")
//' sbs_url <- paste0("https://zenodo.org/records/15875185/files/",
//'   "SBS_demo_signatures.txt")
//' indel_url <- paste0("https://zenodo.org/records/15875185/files/",
//'   "indel_demo_signatures.txt")
//' drivers_url <- paste0("https://zenodo.org/records/15875185/files/",
//'   "driver_mutations_hg19.csv.bz2")
//' passenger_CNAs_url <- paste0("https://zenodo.org/records/15875185/",
//'   "files/passenger_CNAs_hg19.csv.bz2")
//' germline_url <- paste0("https://zenodo.org/records/15875185/files/",
//'   "germline_data_demo.tar.bz2")
//'
//' # build a mutation engine and save the required files into the
//' # directory "Test". The `drivers_url` parameter is optional, but
//' # it is suggested to avoid passenger mutations on driver loci.
//' m_engine <- MutationEngine(directory = "Test",
//'   reference_src = reference_url,
//'   SBS_signatures_src = sbs_url,
//'   indel_signatures_src = indel_url,
//'   drivers_src = drivers_url,
//'   passenger_CNAs_src = passenger_CNAs_url,
//'   germline_src = germline_url)
//'
//' # if the parameters of a mutation engine construction match those of a
//' # previous construction, then the corresponding reference sequence,
//' # the SBS file, and the previously built context index are loaded from
//' # the set-up directory avoiding further computations.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'   drivers_url, passenger_CNAs_url, germline_url)
//'
//' # if the `context_sampling` parameter changes, a new context index is
//' # built, while neither the reference sequence nor the SBS file are
//' # downloaded again.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'   drivers_url, passenger_CNAs_url, germline_url,
//'   context_sampling = 50)
//'
//' # a further construction with the same parameters avoids both
//' # downloads and context index construction.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'   drivers_url, passenger_CNAs_url, germline_url,
//'   context_sampling = 50)
//'
//' m_engine
//'
//' # the parameters `directory`, `reference_src`, `SBS_src`, `drivers_src`,
//' # `passenger_CNAs_src`, and `germline_src` can be avoided by providing
//' # the `setup_code` parameter. The set-up code `demo` is provided among
//' # those available for testing purpose.
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # the `context_sampling` can be used also when a pre-defined set-up
//' # configuration is adopted.
//' m_engine <- MutationEngine(setup_code = "demo", context_sampling = 50)
//'
//' m_engine
//'
//' # remove the "demo" directory
//' unlink("demo", recursive = TRUE)
//'
//' # Some of the pre-defined configurations requires to download the mutational
//' # signatures from the COSMIC site which requires an account (e.g., "GRCh37"
//' # and "GRCh38"). The COSMIC account can be passed to `MutationEngine()` as
//' # follows
//' m_engine <- MutationEngine(setup_code = "demo",
//'   COSMIC_account = list(email = "foo@bar.org",
//'   password = "********"))
//' m_engine
//'
//' # remove the "demo" directory
//' unlink("demo", recursive = TRUE)
//'
//' # In alternative, pre-download the mutational signatures and pass their
//' # paths to `MutationEngine()` as parameters.
//' m_engine <- MutationEngine(setup_code = "demo",
//'   SBS_signatures_src = "Test/SBS_signatures.txt",
//'   indel_signatures_src = "Test/indel_signatures.txt")
//' m_engine
//'
//' # remove the "Test" and "demo" directories
//' unlink("Test", recursive = TRUE)
//' unlink("demo", recursive = TRUE)
    function("MutationEngine", &MutationEngine::build_MutationEngine,
             List::create(_["directory"] = "", _["reference_src"] = "",
                          _["SBS_signatures_src"] = "", _["indel_signatures_src"] = "",
                          _["drivers_src"] = "", _["passenger_CNAs_src"] = "",
                          _["germline_src"] = "", _["setup_code"] = "",
                          _["COSMIC_account"] = R_NilValue, _["germline_subject"] = "",
                          _["context_sampling"] = 30, _["max_motif_size"] = 50,
                          _["max_repetition_storage"] = 500000,
                          _["driver_CNA_min_distance"] = 10000, _["tumour_type"] = "",
                          _["avoid_homozygous_losses"] = true, _["quiet"] = false),
             "Create a MutationEngine");

//' @name get_available_tumours_in
//' @title Getting the tumour types available in a setup
//' @description This method returns the tumour types available for a
//'   set-up code.
//' @usage get_available_tumours_in(setup_code)
//' @param setup_code The set-up code whose available tumour types are
//'   requested.
//' @return A data frame reporting the types available for a set-up code.
//' @seealso [MutationEngine()] to build a mutation engine
//' @export
//' @examples
//' # get the types available for the "demo" set-up code
//' get_available_tumours_in("demo")
    function("get_available_tumours_in", &MutationEngine::get_available_tumour_type,
             List::create(_["setup_code"]),
             "Get the set of tumour types for a set-up code.");

//' @name get_mutation_engine_codes
//' @title Getting the supported setups
//' @description This method returns the supported codes for predefined set-up.
//' @usage get_mutation_engine_codes()
//' @return A data frame reporting the code and a description for each
//'   supported predefined set-up.
//' @seealso [MutationEngine()] to build a mutation engine
//' @export
//' @examples
//' # get the list of supported mutation engine set-up codes
//' get_mutation_engine_codes()
    function("get_mutation_engine_codes", &MutationEngine::get_supported_setups,
             "Get mutation engine supported codes");

//' @name PhylogeneticForest
//' @title The phylogenetic forest of cells in samples
//' @description This class represents the phylogenetic forest of the
//'   cells sampled during the computation.
//' @details The leaves of his forest are the sampled cells.
//'   This class is analogous to the class [SampleForest()],
//'   but each node is labelled with the mutations occurring
//'   for the first time on the cell represented by the node
//'   itself. Moreover each leaf is also associated with the
//'   genome mutations occurring in the corresponding cell.
    class_<PhylogeneticForest>("PhylogeneticForest")

//' @name PhylogeneticForest$get_nodes
//' @title Getting the forest nodes
//' @description This method returns the nodes of the forest.
//' @return A data frame representing, for each node in the
//'   forest, the identified (column `cell_id`), the ancestor
//'   identifier (column `ancestor`), the node's depth (column
//'   `depth`), the name of the sample containing the node
//'   (column `sample`), the mutant (column `mutant`), the
//'   birth time (column `birth_time`), and, whenever the
//'   simulation has epigenetic states, the epigenetic state
//'   (column `epistate`).
//' @seealso [SampleForest$get_nodes()] for usage examples.
        .method("get_nodes", &PhylogeneticForest::get_nodes,
                "Get the nodes of the forest")

//' @name PhylogeneticForest$get_coalescent_cells
//' @title Retrieving the most recent common ancestors
//' @description This method retrieves the most recent common ancestors
//'   of a set of cells.
//' @details If the optional parameter `cell_ids` is
//'   used, this method find the most recent common ancestors of
//'   the cells having an identifier among those in `cell_ids`.
//'   If, otherwise, the optional parameter is not used, this
//'   method find the most recent common ancestors of the forest
//'   leaves.
//' @param cell_ids The list of the identifiers of the cells whose
//'   most recent common ancestors are aimed (optional).
//' @return A data frame representing, for each of the identified
//'   cells, the identified (column `cell_id`), the ancestor
//'   identifier (column `ancestor`), the node's depth (column
//'   `depth`), the name of the sample containing the node
//'   (column `sample`), the mutant (column `mutant`), the
//'   birth time (column `birth_time`), and, whenever the
//'   simulation has epigenetic states, the epigenetic state
//'   (column `epistate`).
//' @seealso [SampleForest$get_coalescent_cells()] for
//'   usage examples
        .method("get_coalescent_cells",
                (List (PhylogeneticForest::*)(const std::list<CLONES::Mutants::CellId> &)
                     const)(&PhylogeneticForest::get_coalescent_cells),
                "Get the most recent common ancestor of some cells")
        .method("get_coalescent_cells",
                (List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_coalescent_cells),
                "Get the most recent common ancestor of all the forest trees")

//' @name PhylogeneticForest$get_subforest_for
//' @title Building sub-forests
//' @description This method builds a sub-forest using as leaves some
//'   of the original samples.
//' @param sample_names The names of the samples whose cells will be used
//'   as leaves of the new forest
//' @return A sample forest built on the samples mentioned in `sample_names`
//' @seealso [SampleForest$get_subforest_for()] for usage examples.
        .method("get_subforest_for", &PhylogeneticForest::get_subforest_for,
                "Get the sub-forest for some of the original samples")

//' @name PhylogeneticForest$partition_samples
//' @title Partitioning forest samples
//' @description This method partitions the samples in a phylogenetic forest.
//' @details This method partitions the samples in a phylogenetic forest
//'   according to a labelling function. It works in-place by altering the
//'   phylogenetic forest from which the method is called.
//' @param labelling_function An R labelling function that maps any sampled
//'   cell to a labelling string.
//' @seealso [SampleForest$get_subforest_for()] for usage examples.
        .method("partition_samples",
                &PhylogeneticForest::partition_samples,
                "Partition each sample in the forest according to a labelling function")

//' @name PhylogeneticForest$get_samples_info
//' @title Retrieving sample information
//' @description This method retrieves information about
//'   the samples whose cells were used as leaves
//'   of the sample forest.
//' @return A data frame containing, for each sample collected during the
//'   simulation, the columns `name`, `time`, `id`, `ymin`, `xmin`,
//'    `ymax`, `ymax`, `xmax`, `tumour_cells`, `tumour_cells_in_bbox`,
//'   `DNA_quantity`, and `equivalent_normal_cells`. The columns `ymin`,
//'   `xmin`, `ymax`, and `xmax` report the boundaries of the sample
//'   bounding box, while `tumour_cells` and `tumour_cells_in_bbox` are the
//'   number of tumour cells in the sample and in the bounding box,
//'   respectively. Finally, `DNA_quantity` contains the overall number of
//'   tumoral bases, i.e., the sum of the lengths of all the alleles of all the
//'   sample tumoral cells, and `equivalent_normal_cells` contains the number
//'   of normal cells that contain as much DNA as the sample tumoral cells.
//'   The `DNA_quantity` is stored as a real number despite being a natural
//'   number as it usually exceeds the largest natural number that can be
//'   natively represented by R.
//' @seealso [SampleForest$get_samples_info()] for usage examples,
//'   [TissueSimulation$get_samples_info()]
        .method("get_samples_info", &PhylogeneticForest::get_samples_info,
                "Get some pieces of information about the samples")

//' @name PhylogeneticForest$get_driver_mutations
//' @title Getting the driver mutations
//' @description This method returns the applied driver mutations.
//' @return A data frame consisting in eight columns `mutant`, `order`, `type`,
//'    `CNA_type`, `chr`, `start`, `end`, `ref`, `alt`, `allele`, `allele_srd`
//'    and `code`. Each row in the data frame reports one driver mutations.
//'    The fields `mutant` and `order` report the name of the mutant and the
//'    application order among the mutant driver mutations, respectively.
//'    The column `type` declares the mutation type and contains `SID`,
//'    `CNA`, or `WGD` when the mutation is an SNV/indel, a CNA, or
//'    a whole genome duplication, respectively. When the mutation is a CNA,
//'    the `CNA_type` can either be `A` (i.e., amplification) or `D`
//'    (i.e., deletion). When the mutation is not a WGD, the fields `chr`,
//'    `start`, and `end` contains the mutation chromosome, the initial and
//'    the final position on the chromosome, respectively. When the mutation
//'    is a SID , the fields `ref` and `alt` contains the mutation reference
//'    genome and alternate sequences, respectively. When the mutation is a
//'    SID or a CNA deletion, the field `allele` stores the allele in which
//'    the mutation was applied. When instead the mutation is a CNA
//'    amplification, the fields `allele` and `src_allele` reports the
//'    identifiers of the new allele and of the original allele, respectively.
//'    In all the remaining cases, the fields contains the value `NA`.
//'    Finally, the column `code` reports the mutation code.
        .method("get_driver_mutations", &PhylogeneticForest::get_driver_mutations,
                "Get the applied driver mutations")

//' @name PhylogeneticForest$get_species_info
//' @title Getting the species and their rates
//' @description This method returns the species and their rates.
//' @details This method returns the species and their rates during
//'   the simulation are returned in a data frame. The column `species`
//'   contains the species names; the columns `time`, `SNV_rate`,
//'   `indel_rate`, and `CNA_rate` store the time from which rates
//'   hold, and the corresponding the SNV, indel, and CNA rates,
//'   respectively.
//' @return A data frame reporting `species`, `time`, `SNV_rate`,
//'   `indel_rate`, and `CNA_rate` for each species.
//' @seealso [MutationEngine$get_species_info()]
        .method("get_species_info", &PhylogeneticForest::get_species_info,
                "Get the recorded species")

//' @name PhylogeneticForest$get_germline_subject
//' @title Getting the germline subject
//' @description This method returns a data frame reporting the germline
//'   subject name (column "sample"), population (column "pop"),
//'   super-population (column "super_pop"), and gender (column "gender").
//' @return The name of the subject whose germline is used.
        .method("get_germline_subject", &PhylogeneticForest::get_germline_subject_df,
                "Get the germline subject")

//' @name PhylogeneticForest$get_sampled_cell_CNAs
//' @title Getting the sampled cells' CNAs
//' @description This method returns the CNAs of the sample cells.
//' @details This method builds a data frame representing all the CNAs
//'   in the cells sampled during the simulation and represented by
//'   the leaves of the phylogenetic forest.
//' @return A data frame reporting `cell_id`, `type` (`"A"` for amplifications
//'   and `"D"` for deletions), `chr`, `begin` (i.e., the first CNA
//'   locus in the chromosome), `end` (i.e., last CNA locus in the chromosome),
//'   `allele`, `src allele` (the allele origin for amplifications, `NA` for
//'   deletions), and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"`
//'   or `"pre-neoplastic"`).
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # build a tissue simulation with epigenetic states "E1" and "E2"
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place a cell of species "A[E1]" in position (500,500)
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation until "A[E2]" accounts for less than 1000 cells
//' sim$run_up_to_size("A[E2]", 1000)
//'
//' # sample the tissue
//' sim$sample_cells("Sample_1", c(475, 475), c(525, 525))
//'
//' # build the sample forest
//' sample_forest <- sim$get_sample_forest()
//'
//' # initialize a mutation engine with the "demo" setup
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add the genomic characterisation for the mutant "A"
//' m_engine$add_mutant("A",
//'                     list("E1" = c(SNV = 1e-7, indel = 1e-8),
//'                          "E2" = c(SNV = 3e-7, CNA = 1e-11)),
//'                     list(SNV("22", 10510210, "C", allele = 1),
//'                          CNA("D", "22", 5010000, 200000,
//'                              allele = 1)))
//'
//' # add the exposure
//' m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))
//'
//' # build the phylogenetic forest
//' phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
//'
//' # get the CNAs in the sampled cells
//' CNAs <- phylo_forest$get_sampled_cell_CNAs()
//'
//' # show the first CNAs in the sampled cells
//' CNAs %>% head()
//' @seealso [PhylogeneticForest$get_sampled_cell_mutations()]
        .method("get_sampled_cell_CNAs", &PhylogeneticForest::get_sampled_cell_CNAs,
                "Get the CNAs of all sampled cells")

//' @name PhylogeneticForest$get_sampled_cell_mutations
//' @title Getting the sampled cells' mutations
//' @description This method returns the mutations of the sample cells.
//' @details This method builds a data frame representing all the SNV
//'   and the indel mutations in the cells sampled during the simulation
//'   and represented by the leaves of the phylogenetic forest.
//'   The data frame also reports the allele in which the mutations occur to
//'   support double occurrences due to CNAs.
//' @param with_germline A Boolean flag to report germline mutations too
//'   (default: `FALSE`).
//' @return A data frame reporting `cell_id`, `chr`, (i.e., the mutation
//'   chromosome), `from` (i.e., position in the chromosome), `allele`
//'   (in which the mutation occurs), `ref`, `alt`, `type` (i.e., either
//'   `"SNV"` or `"indel"`), `cause`, and `class` (i.e., `"driver"`,
//'   `"passenger"`, `"germinal"` or `"pre-neoplastic"`) for each mutation
//'   in the sampled cell genomes.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # build a tissue simulation with epigenetic states "E1" and "E2"
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place a cell of species "A[E1]" in position (500,500)
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation until "A[E2]" accounts for less than 1000 cells
//' sim$run_up_to_size("A[E2]", 1000)
//'
//' # sample the tissue
//' sim$sample_cells("Sample_1", c(475, 475), c(525, 525))
//'
//' # build the sample forest
//' sample_forest <- sim$get_sample_forest()
//'
//' # initialize a mutation engine with the "demo" setup
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add the genomic characterisation for the mutant "A"
//' m_engine$add_mutant("A",
//'                     list("E1" = c(SNV = 1e-7, indel = 1e-8),
//'                          "E2" = c(SNV = 3e-7, CNA = 1e-11)),
//'                     list(SNV("22", 10510210, "C", allele = 1),
//'                          CNA("D", "22", 5010000, 200000,
//'                              allele = 1)))
//'
//' # add the exposure
//' m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))
//'
//' # build the phylogenetic forest
//' phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
//'
//' # get the SIDs in the sampled cells
//' SIDs <- phylo_forest$get_sampled_cell_mutations()
//'
//' # show the first SIDs in the sampled cells
//' SIDs %>% head()
//' @seealso [PhylogeneticForest$get_sampled_cell_CNAs()]
        .method("get_sampled_cell_mutations",
                (Rcpp::DataFrame (PhylogeneticForest::*)() const)(
                    &PhylogeneticForest::get_sampled_cell_SIDs),
                "Get the SNVs and indels of all the sampled cells")
        .method("get_sampled_cell_mutations",
                (Rcpp::DataFrame (PhylogeneticForest::*)(const bool) const)(
                    &PhylogeneticForest::get_sampled_cell_SIDs),
                "Get the SNVs and indels of all the sampled cells")

//' @name PhylogeneticForest$get_node
//' @title Getting a node of the forest
//' @description This method returns the node of the forest
//' @details This method returns the node of the forest whose
//'    corresponding cell has a specified identifier.
//' @param cell_id The identifier of the cell whose node is aimed.
//' @return The <code>[PhylogeneticForestNode]</code> object
//'   associated to the cell whose identifier is `cell_id`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # build a tissue simulation with epigenetic states "E1" and "E2"
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add mutant "A" and set its species rates
//' sim$add_mutant("A",
//'                list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
//'                     E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))
//'
//' # place a cell of species "A[E1]" in position (500,500)
//' sim$place_cell("A[E1]", 500, 500)
//'
//' # run the simulation until "A[E2]" accounts for less than 1000 cells
//' sim$run_up_to_size("A[E2]", 1000)
//'
//' # sample the tissue
//' sim$sample_cells("Sample_1", c(475, 475), c(525, 525))
//'
//' # build the sample forest
//' sample_forest <- sim$get_sample_forest()
//'
//' # initialize a mutation engine with the "demo" setup
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add the genomic characterisation for the mutant "A"
//' m_engine$add_mutant("A",
//'                     list("E1" = c(SNV = 1e-7, indel = 1e-8),
//'                          "E2" = c(SNV = 3e-7, CNA = 1e-11)),
//'                     list(SNV("22", 10510210, "C", allele = 1),
//'                          CNA("D", "22", 5010000, 200000,
//'                              allele = 1)))
//'
//' # add the exposure
//' m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))
//'
//' # build the phylogenetic forest
//' phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
//'
//' # get a node corresponding to a non-sampled cell
//' node <- phylo_forest$get_nodes() %>%
//'   dplyr::filter(is.na(sample)) %>%
//'   dplyr::sample_n(1)
//'
//' # gets the SIDs in it
//' phylo_forest$get_node(node$cell_id)
//' @seealso [PhylogeneticForest$get_sampled_cell_mutations()],
//'   [PhylogeneticForest$get_sampled_cell_CNAs()]
        .method("get_node", &PhylogeneticForest::get_node,
                "Get node corresponding to a cell")

//' @name PhylogeneticForest$get_germline_mutations
//' @title Getting the germinal mutations
//' @description This method returns the forest SNVs and indels.
//' @details Its builds a data frame representing all the germinal
//'   SNVs and indels of the cells represented in the phylogenetic forest.
//'   The data frame also reports the allele in which the mutations occur to
//'   support double occurrences due to CNAs.
//' @return A data frame reporting `chr`, `from` (i.e., the position in
//'   the chromosome), `allele` (in which the mutation occurs), `ref`,
//'   `alt`, `cause`, `type` (i.e., either `"SNV"` or `"indel"`) and
//'   `class` (i.e., `"germinal"`).
//' @seealso [`vignette("mutations")`] for usage examples.
        .method("get_germline_mutations", &PhylogeneticForest::get_germline_SIDs,
                "Get the germinal SNVs and indels")

//' @name PhylogeneticForest$get_absolute_chromosome_positions
//' @title Getting the absolute chromosome positions
//' @description This method returns the absolute chromosome positions.
//' @details Its builds a data frame reporting the name, the length, and the
//'   initial and final absolute positions of each chromosome in the
//'   reference genome.
//' @return A data frame reporting the name (column `chr`), the length
//'   (column `length`), the initial absolute position (column `from`),
//'   and the final absolute position (column `to`) of each chromosome.
        .method("get_absolute_chromosome_positions",
                &PhylogeneticForest::get_absolute_chromosome_positions,
                "Get the absolute chromosome positions")

//' @name PhylogeneticForest$get_sticks
//' @title Computing the forest sticks
//' @description This method computes the sticks of the forest.
//' @details A _crucial node_ of a forest is a root of the forest, a node
//'   whose parent belongs to a different species, or the most recent common
//'   ancestor of two crucial nodes.
//'
//'   A _stick_ is a path of the forest in which the only crucial nodes are
//'   the first and the last one.
//'
//'   This method returns the list of the forest sticks. Each stick is
//'   represented by the sequence of cell identifiers labelling the nodes in
//'   the stick.
//' @param birth_threshold The maximum birth time for the cells associated to
//'   the returned sticks (optional).
//' @return The list of the forest sticks whose associated cells have
//'   birth time smaller than or equal to `birth_threshold`. Each stick is
//'   represented as the list of cell identifiers labelling the nodes in the
//'   stick from the higher to the deeper in the forest.
//' @seealso [SampleForest$get_sticks()] for usage examples.
        .method("get_sticks",
                (std::list<std::list<CLONES::Mutants::CellId>> (PhylogeneticForest::*)(
                    const double) const)(&PhylogeneticForest::get_sticks),
                "Get the forest sticks")
        .method("get_sticks",
                (std::list<std::list<CLONES::Mutants::CellId>> (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_sticks),
                "Get the forest sticks")

//' @name PhylogeneticForest$get_exposures
//' @title Getting the timed exposure data frame
//' @description This method returns a data frame representing the exposure
//'   evolution over time.
//' @return A data frame reporting `time`, `signature`, `exposure` and,
//'   `type`.
//' @seealso [`vignette("mutations")`] for usage examples.
        .method("get_exposures", &PhylogeneticForest::get_timed_exposures,
                "Get the timed exposure data frame")

//' @name PhylogeneticForest$get_bulk_allelic_fragmentation
//' @title Getting the bulk allelic fragmentation data frame
//' @description This method returns a data frame representing the bulk allelic
//'   fragmentation of the genome.
//' @param sample_name The name of the sample whose bulk allelic fragmentation
//'   is aimed (optional).
//' @return A data frame reporting, for each genomic fragment and for all
//'   the allelic type on the genomic fragment, the chromosome (`chr`),
//'   the first base position (`begin`), the last base position (`end`),
//'   the number of copy of the major and minor alleles (`major` and
//'   `minor`, respectively), and the ratio between the number of cells
//'   exhibiting this allelic type and the total number of cells in the
//'   sample.
//' @seealso [`vignette("mutations")`] for usage examples.
        .method("get_bulk_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)(const std::string &)
                     const)(&PhylogeneticForest::get_bulk_allelic_fragmentation),
                "Get the bulk allelic fragmentation data frame")
        .method("get_bulk_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_bulk_allelic_fragmentation),
                "Get the bulk allelic fragmentation data frame")

//' @name PhylogeneticForest$get_cell_allelic_fragmentation
//' @title Getting the cell allelic fragmentation data frame
//' @description This method returns a data frame representing the allelic
//'   fragmentation of each sampled cell.
//' @return A data frame reporting, for each cell, for each genomic fragment,
//'   and for all the allelic type on the genomic fragment, the cell
//'   identifier (`cell_id`), the chromosome (`chr`), the first base
//'   position (`begin`), the last base position (`end`), and the number
//'   of copy of the major and minor alleles (`major` and `minor`,
//'   respectively).
//' @seealso [`vignette("mutations")`] for usage examples.
        .method("get_cell_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_cell_allelic_fragmentation),
                "Get the cell allelic fragmentation data frame")

//' @name PhylogeneticForest$get_first_occurrences
//' @title Getting the cell in which a mutation emerged
//' @description This method returns the identifier of the cell in which a
//'   mutation occurs for the first time. `$1`
//' @param mutation A mutation being a SNV, a indel, or a CNA.
//' @return The identifier of the cell in which a mutation
//'   occurs for the first time.
//' @seealso [`vignette("mutations")`] for usage examples.
        .method("get_first_occurrences",
                (Rcpp::List (PhylogeneticForest::*)(const SEXP &)
                     const)(&PhylogeneticForest::get_first_occurrence),
                "Get the identifier of the cell in which the mutation occurs for the "
                "first time")

//' @name PhylogeneticForest$get_reference_path
//' @title Getting the reference genome path
//' @description This method returns the reference genome path.
//' @return The reference genome path.
//' @seealso [PhylogeneticForest$set_reference_path()]
        .method("get_reference_path",
                &PhylogeneticForest::get_reference_path,
                "Get the reference genome path")

//' @name PhylogeneticForest$get_mutation_statistics
//' @title Getting the statistics about mutations on each node
//' @description This method returns a dataframe reporting the statistics about
//'   mutations on each node.
//' @return A dataframe consisting of five columns `cell_id`, `new_SIDs`,
//'   `new_CNAs`, `total_SIDs`, and `total_CNAs`. Each row represents a
//'   node in the phylogenetic forest and reports the identifier of the
//'   corresponding cell and contains the number of mutations (`new_SIDs`) and
//'   CNAs (`new_CNAs`) appearing for the first time on the cell. Moreover, it
//'   show the total number of mutations and CNAs on the cell (`total_SIDs`
//'   and `total_CNAs`, respectively).
        .method("get_mutation_statistics",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_mutation_statistics),
                "Get the statistics about node mutations")

//' @name PhylogeneticForest$set_reference_path
//' @title Setting the reference genome path
//' @description This method returns the reference genome path.
//' @return The reference genome path.
//' @seealso [PhylogeneticForest$get_reference_path()]
        .method("set_reference_path",
                &PhylogeneticForest::set_reference_path,
                "Set the reference genome path")

//' @name PhylogeneticForest$save
//' @title Saving a phylogenetic forest
//' @description This method saves a phylogenetic forest in a file.
//' @param filename The path of the file in which the phylogenetic.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: FALSE).
//'   forest must be saved.
        .method("save",
                (void (PhylogeneticForest::*)(const std::string &, const bool) const)
                    &PhylogeneticForest::save,
                "Save a phylogenetic forest")
        .method("save",
                (void (PhylogeneticForest::*)(const std::string &) const)
                    &PhylogeneticForest::save,
                "Save a phylogenetic forest")

        // show
        .method("show", &PhylogeneticForest::show, "Describe the PhylogeneticForest");

//' @name PhylogeneticForestNode
//' @title The node of a phylogenetic forest
//' @description This class represents the nodes of a phylogenetic forest. It does not
//'   have a user constructor and its objects are produced by ProCESS and passed to
//'   the labelling function of [get_node_tour()].
//' @field \code{cell_id} The identifier of the associated cell.
//' @field \code{parent} The node's parent.
//' @field \code{children} A list of the node's children.
//' @field \code{is_root} A flag that is set to TRUE if and only if the node is a root.
//' @field \code{is_leaf} A flag that is set to TRUE if and only if the node is a leaf.
//' @field \code{sample_name} The name of the sample that collected the associated cell.
//' @field \code{birth_time} The birth time of the cell associated to the node.
//' @field \code{death_time} The death time of the cell associated to the node.
//' @field \code{life_span} The life span of the cell associated to the node.
//' @field \code{species_id} The identifier of the associated cell's species.
//' @field \code{species_name} The name of the associated cell's species.
//' @field \code{epistate_name} The name of the associated cell's epigenetic state.
//' @field \code{mutant_id} The identifier of the associated cell's mutant.
//' @field \code{mutant_name} The name of the associated cell's mutant.
//' @field \code{arising_mutations} The mutations arising for the first time in the
//'   associated cell.
//' @seealso [get_node_tour()], <code>[PhylogeneticForestNodeTour]</code>,
//'   <code>[SampleForestNode]</code>, [`vignette("node_labelling")`]
   class_<PhylogeneticForestNode>("PhylogeneticForestNode")
        REGISTER_NODE_COMMON_FIELDS(PhylogeneticForestNode)
        .property("arising_mutations", &PhylogeneticForestNode::arising_mutations,
            "The mutations occurring for the first time in the associated cell (Read-only)")
        .method("get_genome", &PhylogeneticForestNode::get_genome,
            "Get the genome of the associated cell");

//' @name PhylogeneticForestNodeTour
//' @title An iterator class over phylogenetic forest nodes
//' @description This class represents iterators over phylogenetic forest nodes.
//'   The objects of this class are built by [get_node_tour()].
//' @field \code{node} An object of the class <code>[PhylogeneticForestNode]</code>
//'   representing the node pointed by the iterator.
//' @field \code{label} (OPTIONAL) The label of the of the node pointed by the
//'   iterator. The presence of this field depends on the [get_node_tour()]'s
//'   parameters used to create the tour object.
//' @field \code{genome} (OPTIONAL) An object of the class
//'   <code>[GenomeMutations]</code> that represent the genome of the node
//'   pointed by the iterator.
//' @field \code{step} A method that moves to the next node in the tour.
//' @field \code{done} A Boolean flag that is set to TRUE only when the tour ended.
//' @seealso [get_node_tour()], <code>[PhylogeneticForestNode]</code>,
//'   <code>[SampleForestNodeTour]</code>, <code>[GenomeMutations]</code>,
//'   [`vignette("node_labelling")`]
    class_<PhylogeneticForestNodeTour>("PhylogeneticForestNodeTour")
        REGISTER_TOUR_COMMON_FIELDS(PhylogeneticForestNodeTour)
        .property("genome", &PhylogeneticForestNodeTour::get_genome,
            "The genome of the node pointed by the iterator");

//' @name GenomeFragment
//' @title Representing a genome fragment
//' @description This class represents genome fragment.
//'   The objects of this class are built by <code>[GenomeMutations]</code>.
//' @seealso <code>[GenomeMutations]</code>, [get_node_tour()],
//'   [get_genome_tour()], [`vignette("node_labelling")`]
    class_<GenomeFragment>("GenomeFragment")

//' @name GenomeFragment$get_CIGAR
//' @title Getting the fragment CIGAR
//' @description This method returns the CIGAR code of the fragment with respect
//'   to the reference genome.
//' @return The CIGAR code of the fragment with respect to the reference genome.
        .method("get_CIGAR", &GenomeFragment::get_CIGAR,
                "Get the fragment CIGAR")

//' @name GenomeFragment$sequence
//' @title The fragment sequence
//' @description This property is the fragment sequence.
//' @return The fragment sequence.
        .property("sequence", &GenomeFragment::get_sequence,
                  "The fragment sequence")

//' @name GenomeFragment$get_mutations
//' @title Getting the mutations laying on a fragment
//' @description This method returns the data frame of the mutations
//'   laying on the fragment.
//' @return A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`, `alt`,
//'   `cause`, and `nature`. Each row represent a SID.
        .method("get_mutations", &GenomeFragment::get_mutations,
                "Get the mutations laying on the fragment")

//' @name GenomeFragment$get_covered_reference_region
//' @title Getting the reference region covered by the fragment
//' @description This method returns the reference region covered by the fragment.
//' @return A named list representing the reference region covered by the fragment.
//'   The names of the list are `chr`, `from`, and `size`. They stores the name
//'   of the chromosome from which the fragment comes, the position in the
//'   chromosome of the fragment first base, and the size of the covered region,
//'   respectively.
        .method("get_covered_reference_region",
                &GenomeFragment::get_covered_reference_region,
                "Get the covered region")

        .method("show", &GenomeFragment::show);

//' @name GenomeMutations
//' @title Representing cell genome
//' @description This class represents genome mutations of phylogenetic forest's cells.
//'   The objects of this class are built by <code>[PhylogeneticForestNodeTour]</code>.
//' @seealso [get_node_tour()], [get_genome_tour()],
//'   <code>[PhylogeneticForestNodeTour]</code>, [`vignette("node_labelling")`]
    class_<GenomeMutations>("GenomeMutations")

//' @name GenomeMutations$get_mutations
//' @title Getting genome mutations
//' @description This method returns a data frame representing the SIDs in the genome.
//' @param with_germline A Boolean flag to add germline mutations (default: `TRUE`).
//' @return A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`, `alt`,
//'   `cause`, and `nature`. Each row represent a SID.
        .method("get_mutations",
            (Rcpp::DataFrame (GenomeMutations::*)() const)(&GenomeMutations::get_mutations),
            "Get the genome mutations")
        .method("get_mutations",
            (Rcpp::DataFrame (GenomeMutations::*)(const bool) const)(&GenomeMutations::get_mutations),
            "Get the genome mutations")

//' @name GenomeMutations$get_CNAs
//' @title Getting genome mutations
//' @description This method returns a data frame representing the CNAs in the genome.
//' @return A data frame consisting of 8 columns: `chr`, `begin`, `end`, `type`,
//'   `allele`, `src.allele`, `cause`, and `nature`. Each row represent a SID.
        .method("get_CNAs",
            (Rcpp::DataFrame (GenomeMutations::*)() const)(&GenomeMutations::get_CNAs),
            "Get the genome CNAs")

//' @name GenomeMutations$get_allele_fragments
//' @title Getting genome allele fragments
//' @description This method returns a data frame representing the allele fragments
//'   in the genome.
//' @return A data frame consisting of 5 columns: `chr`, `allele`, `src_allele`,
//'   `from`, and `size`. Each row represent a allele fragment. The columns `chr`,
//'   and `allele` represent the fragment's chromosome and allele, respectively.
//'   The column `allele_src` stores the allele from which the allele of the
//'   fragment is derived. The columns `from` and `size` maintain the first
//'   position of the fragment in the wild-type allele and its size.
        .method("get_allele_fragments",
            (Rcpp::DataFrame (GenomeMutations::*)() const)(&GenomeMutations::get_allele_fragments),
            "Get the genome allele fragments")

//' @name GenomeMutations$get_region_aligned_to_ref
//' @title Getting information about the fragment aligning to the reference
//' @description This method returns information about the fragment of the
//'   current genome that align with the provided region in the reference.
//' @param chromosome_name The name of the chromosome of the reference region.
//' @param allele The allele of the reference region.
//' @param from The first position in the chromosome of the reference region.
//' @param size The size of the reference region.
//' @return A named list whose names are `chr`, `allele`, `from`, and `length`
//'   representing the fragment in the current genome that aligns on
//'   the region in chromosome `chromosome_name`, allele `allele` from
//'   position `from` and whose size is `size`.
        .method("get_region_aligned_to_ref",
                &GenomeMutations::region_aligning_on_reference,
            "Get information about the fragment aligning to the specified "
            "reference region")

//' @name GenomeMutations$get_fragment
//' @title Getting genome fragment
//' @description This method returns a fragment of the genome.
//' @param chromosome_name The name of the chromosome of the aimed fragment.
//' @param allele The allele of the aimed fragment.
//' @param from The first position in the chromosome of the aimed fragment.
//' @param size The size of the aimed fragment.
//' @param reference_fragment A reference fragment (optional).
//' @param fragment_offset The offset of the reference fragment with respect
//'   to the chromosome first position (optional).
//' @return The genome fragment matching the specifications.
//' @seealso <code>[GenomeFragment]</code>
//'   <code>[GenomeMutations$get_region_aligned_to_ref()]</code>
        .method("get_fragment",
                (GenomeFragment (GenomeMutations::*)(const std::string& chromosome_name,
                                                     const size_t& allele_id,
                                                     const size_t& from,
                                                     const size_t& size) const)
                (&GenomeMutations::get_fragment),
            "Get a genome fragment")
        .method("get_fragment",
                (GenomeFragment (GenomeMutations::*)(const std::string& chromosome_name,
                                                     const size_t& allele_id,
                                                     const size_t& from,
                                                     const size_t& size,
                                                     const std::string& reference_fragment,
                                                     const size_t& fragment_offset) const)
                (&GenomeMutations::get_fragment),
            "Get a genome fragment")

        .method("show", &GenomeMutations::show);

//' @name load_phylogenetic_forest
//' @title Loading a phylogenetic forest
//' @description This method loads a phylogenetic forest from a file.
//' @usage load_phylogenetic_forest(filename)
//' @param filename The path of the file from which the phylogenetic
//'   forest must be load.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @return The load phylogenetic forest
    function("load_phylogenetic_forest",
             (PhylogeneticForest (*)(const std::string &,
                                     const bool))&PhylogeneticForest::load,
             List::create(_["filename"] = "", _["quiet"] = false),
             "Load a phylogenetic forest");
}
