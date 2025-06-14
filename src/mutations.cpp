/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2025 Alberto Casagrande <alberto.casagrande@uniud.it>
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

#include "phylogenetic_forest.hpp"
#include "mutation_engine.hpp"
#include "sid.hpp"
#include "cna.hpp"
#include "wg_doubling.hpp"

using namespace Rcpp;

RCPP_MODULE(Mutations){

//' @name Mutation
//' @title A single nucleotide variation
  class_<SIDMut>("Mutation")
    .constructor()

//' @name Mutation$get_chromosome
//' @title Getting the mutation chromosome
//' @description This method returns the identifier of the chromosome
//'   where the mutation occurs.
//' @return The chromosome in which the mutation occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the chromosome in which `snv` occurs (i.e., "X")
//' snv$get_chromosome()
    .method("get_chromosome",&SIDMut::get_chromosome,
            "Get the chromosome of the mutation")

//' @name Mutation$get_position_in_chromosome
//' @title Getting the mutation chromosome position
//' @description This method returns the position in the chromosome where
//'   the mutation occurs.
//' @return The position in chromosome where the mutation occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the position in chromosome where `snv` occurs (i.e., 20002)
//' snv$get_position_in_chromosome()
    .method("get_position_in_chromosome",&SIDMut::get_position_in_chromosome,
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
    .method("get_ref", &SIDMut::get_ref,
            "Get the mutation reference sequence")

//' @name Mutation$get_alt
//' @title Getting the mutation altered sequence
//' @description This method returns the sequence after the mutation occurs.
//' @return The sequence after the mutation occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the sequence after `snv` occurs (i.e., "T")
//' snv$get_alt()
    .method("get_alt", &SIDMut::get_alt,
            "Get the mutation altered sequence")

//' @name Mutation$get_cause
//' @title Getting the mutation cause
//' @description This method returns the mutation cause.
//' @details Evey mutation is associated to a cause depending on whether
//'   it is part of a genomic characterization of a mutant or it is caused
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
//' @title Getting the mutation dataframe
//' @description This method builds a dataframe representing the mutation.
//' @details The dataframe has the columns `chr`, `chr_pos`, `ref`, `alt`,
//'   `type` (i.e., "`SNV`" and "`indel`"), and `cause`.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' snv$get_dataframe()
    .method("get_dataframe", &SIDMut::get_dataframe,
            "Get a dataframe representing the mutation")
    .method("show", &SIDMut::show);

//' @name SNV
//' @title Creating a SNV
//' @description This function creates SNVs.
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
  function("SNV",  &SIDMut::build_SNV,
           List::create(_["chr"], _["chr_pos"], _["alt"],
                        _["ref"] = "?", _["allele"] = R_NilValue,
                        _["cause"] = ""),
           "Create a single nucleotide variation (SNV)");

//' @name Mutation
//' @title Creating a SNV or a indel
//' @description This function creates SNVs and indels.
//' @details It generalizes the function `SNV()` by building SNVs and
//'   indels. However, it requires the reference sequence specification
//'   whereas `SNV()` can deduce it from the reference sequence itself.
//'
//'   Another difference with respect to `SNV()` is the `ref`-`alt`
//'   parameter order: the `alt` parameter comes before the optional
//'   `ref` parameter in `SNV()`; `Mutation()` adopts the reverse order.
//' @param chr The name of the chromosome in which the indel occurs.
//' @param chr_pos The position in the chromosome where the indel occurs.
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
  function("Mutation",  &SIDMut::build_SID,
           List::create(_["chr"], _["chr_pos"], _["ref"],
                        _["alt"], _["allele"] = R_NilValue,
                        _["cause"] = ""),
           "Create either an SNV or a indel");

//' @name CNA
//' @title Creating a CNA
//' @description This function creates a CNA.
//' @param type The CNA type: either `"A"` or `"D"` for amplification and
//'   deletion, respectively.
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param chr_pos The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the CNA occurs. (optional)
//' @param src_allele The allele from which the region is amplified. (optional,
//'   for amplification only)
//' @seealso `Amplification` to build an amplification; `Deletion` to build a
//'   deletion.
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
           List::create(_["type"], _["chr"], _["chr_pos"], _["len"],
                        _["allele"] = R_NilValue, _["src_allele"] = R_NilValue),
           "Create a copy number alteration (CNA)");

//' @name Amplification
//' @title Creating a CNA amplification
//' @description This function creates a CNA amplification.
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param chr_pos The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the amplification is placed. (optional)
//' @param src_allele The allele from which the region is amplified. (optional)
//' @seealso `Deletion` to build a deletion; `CNA` to build both amplifications
//'   and deletions.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100)
//'
//' cna
  function("Amplification", &CNA::build_amplification,
           List::create(_["chr"], _["chr_pos"], _["len"],
                        _["allele"] = R_NilValue, _["src_allele"] = R_NilValue),
           "Create a CNA amplification");

//' @name Deletion
//' @title Creating a CNA deletion
//' @description This function creates a CNA deletion.
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param chr_pos The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the deletion occurs. (optional)
//' @seealso `Amplification` to build an amplification; `CNA` to build
//'   both amplifications and deletions.
//' @examples
//' # create a deletion CNA
//' cna <- Deletion("Y", 40020, 200)
//'
//' cna
  function("Deletion", &CNA::build_deletion,
           List::create(_["chr"], _["chr_pos"], _["len"],
                        _["allele"] = R_NilValue),
           "Create a CNA deletion");

//' @name CNA
//' @title A copy number alteration
  class_<CNA>("CNA")
    .constructor()

//' @name CNA$get_chromosome
//' @title Getting the CNA chromosome
//' @description This method returns the identifier of the chromosome
//'   where the CNA occurs.
//' @return The identifier of the chromosome in which the CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- CNA("A", "X", 20002, 100)
//'
//' # get the chromosome in which `cna` occurs (i.e., "X")
//' cna$get_chromosome()
    .method("get_chromosome",&CNA::get_chromosome, "Get the chromosome of the CNA")

//' @name CNA$get_position_in_chromosome
//' @title Getting the CNA chromosome position
//' @description This method returns the position in chromosome
//'   where the CNA occurs.
//' @return The position in chromosome where the CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the position in chromosome where `cna` occurs (i.e., 20002)
//' cna$get_position_in_chromosome()
    .method("get_position_in_chromosome",&CNA::get_position_in_chromosome,
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
    .method("get_length",&CNA::get_length, "Get the CNA length")

//' @name CNA$get_allele
//' @title Getting the CNA allele
//' @description This method returns the identifier of the allele in
//'    which the CNA is occurs.
//' @details If the CNA is an amplification corresponds to the new
//'    allele identifier. If, instead, the CNA is a deletion is the
//'    identifier of the allele on which the deletion occurs.
//' @return The allele in which CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the allele in which `cna` occurs (i.e., 1)
//' cna$get_allele()
    .method("get_allele",&CNA::get_allele, "Get the alteration allele")

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
    .method("get_src_allele",&CNA::get_src_allele, "Get the source allele (for amplifications)")

//' @name CNA$get_dataframe
//' @title Getting the CNA dataframe
//' @description This method builds a dataframe representing the CNA.
//' @details The dataframe contains the  columns "`chr`", "`chr_pos`",
//'   "`length`", "`alt_base`", "`allele`"", "`src_allele`", and "`type`".
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
    .method("get_dataframe",&CNA::get_dataframe, "Get a dataframe representing the CNA")
    .method("show",&CNA::show);

//' @name WholeGenomeDoubling
//' @title Whole genome doubling events
//' @description A whole genome doubling event (WGD)
//'   produces the simultaneous duplication of all the
//'   chromosome allele in a genome.
  class_<WholeGenomeDoubling>("WholeGenomeDoubling")
    .constructor()
    .method("show",&WholeGenomeDoubling::show);

//' @name MutationEngine
//' @title Generating phylogenetic forests
//' @description A mutation engine can label every node of a descendants
//'   forest by mutations and produce a consistent phylogenetic forest.
//' @details The mutations are randomly generated according to three
//'   factors:
//'     - the mutational rates of the species involved in the descendants
//'       forest
//'     - the genotypical characterisation of the mutants involved in the
//'       descendants forest, i.e., the somatic mutations characterising
//'       the mutant genotypes
//'     - the SBS and ID signatures active along the species simulation
//'
//'   These data are provided to a mutation engine by using the methods
//'   [MutationEngine$add_exposure()] and [MutationEngine$add_exposure()]
//'   These data are provided by means of the [MutationEngine$add_mutant()].
//'
//'   The initialisation of a `MutationEngine` object requires a reference
//'   sequence and the SBS and ID mutational signatures. An SBS index and
//'   a ID index of the reference sequence are then automatically built.
//'   Thess processes may take time depending on the size of the reference
//'   sequence. Because of this, the downloaded files together with the
//'   context index are saved in a directory on the disk and they are
//'   available for successive `MutationEngine` constructions.
//'
  class_<MutationEngine>("MutationEngine")

//' @name MutationEngine$infinite_sites_model
//' @title Switching on and off the infinite sites model.
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
//'                         ID20 = 0.5))
//'
//' # add a default set of coefficients that will be used from simulated
//' # time 3.2 up to the end of the simulation.
//' m_engine$add_exposure(3.2, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))
//'
//' m_engine
    .method("add_exposure", (void (MutationEngine::*)(const List&))(
              &MutationEngine::add_exposure), "Add an exposure")
    .method("add_exposure", (void (MutationEngine::*)(const double&,
                                                      const List&))(
              &MutationEngine::add_exposure), "Add an exposure")

//' @name MutationEngine$add_mutant
//' @title Adding a mutant specification
//' @description This method adds a mutant specification to the mutation engine.
//' @details The users must use it to specify the name and the genomic
//'   characterization (i.e., SNVs, indels, CNAs, and whole genome doubling
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
//' @examples
//' # create a demonstrative mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # add the mutant "A" characterized by one driver SNV on chromosome 22, an
//' # indel on the same chromosome, a whole genome doubling event, and finally
//' # two CNAs: an amplification and a deletion. The mutant has two epigenetic
//' # states and its species "A+" and "A-" have passenger SNV rates 1e-9 and
//' # 3e-8, respectively, and passenger CNA rates 0 and 1e-11, respectively.
//' m_engine$add_mutant("A", list("+" = c(SNV = 1e-9, indel = 1e-10),
//'                               "-" = c(SNV = 3e-8, CNA = 1e-11)),
//'                     drivers = list(SNV("22", 23657587, "C"),
//'                                    Mutation("22", 15220157, "GTTTTTTTT",
//'                                             "G"),
//'                                    WGD,
//'                                    CNA(type = "A", chr = "22",
//'                                        chr_pos = 10303470,
//'                                        len = 200000),
//'                                    CNA("D", "22", 5010000, 200000)))
//'
//' m_engine
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& passenger_rates))(
                                                        &MutationEngine::add_mutant),
            "Add mutant")
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& passenger_rates,
                                                    const Rcpp::List&))(
                                                        &MutationEngine::add_mutant),
            "Add mutant")

//' @name MutationEngine$place_mutations
//' @title Placing the mutations
//' @description This methods places mutations on a samples forest.
//' @details Each node of a samples forest is labelled by the
//'   mutations occurring in the cell represented by the node itself
//'   and produces a phylogenetic forest.
//' @param samples_forest A samples forest.
//' @param num_of_preneoplatic_SNVs The number of pre-neoplastic SNVs.
//' @param preneoplatic_SNV_signature_name The name of the SNV signature
//'   for the preneoplastic SNV generation (optional).
//' @param num_of_preneoplatic_indels The number of pre-neoplastic indels.
//' @param preneoplatic_indel_signature_name The name of the indel signature
//'   for the preneoplastic indel generation.
//' @param seed The seed for random number generator (optional).
//' @return A phylogenetic forest whose structure corresponds to
//'   `samples_forest`.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant("A", c(SNV = 0.2), 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner = c(450, 475),
//'                        upper_corner = c(500, 550))
//'
//' # build the samples forest
//' samples_forest <- sim$get_samples_forest()
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
//' # place the mutations on the samples forest assuming 1000 pre-neoplastic
//' # SNVs and 500 indels
//' phylogenetic_forest <- m_engine$place_mutations(samples_forest, 1000, 500)
//'
//' phylogenetic_forest
    .method("place_mutations",
            (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                    const size_t& num_of_preneoplatic_SNVs,
                                                    const size_t& num_of_preneoplatic_indels))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("place_mutations",
            (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                    const size_t& num_of_preneoplatic_SNVs,
                                                    const size_t& num_of_preneoplatic_indels,
                                                    const SEXP& seed))(&MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("place_mutations",
            (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                    const size_t& num_of_preneoplatic_SNVs,
                                                    const std::string& preneoplatic_SNV_signature_name,
                                                    const size_t& num_of_preneoplatic_indels,
                                                    const std::string& preneoplatic_indel_signature_name))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("place_mutations",
            (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                    const size_t& num_of_preneoplatic_SNVs,
                                                    const std::string& preneoplatic_SNV_signature_name,
                                                    const size_t& num_of_preneoplatic_indels,
                                                    const std::string& preneoplatic_indel_signature_name,
                                                    const SEXP& seed))(&MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")

//' @name MutationEngine$get_genome_info
//' @title Getting the genome information
//' @description This method returns information about the genome.
//' @details This method returns a data frame reporting the name
//'    (column "`name`"), the size (column "`size`"), and the number
//'    of alleles (column "`num_of_alleles`") of each chromosome.
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
//'   dataframe in which the column `sample` reports the
//'   subject name, the columns `pop` and `super_pop` contain the
//'   subject population and super population, respectively, and
//'   the column `gender` declares the subject gender.
//' @return A dataframe the active germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to get the
//'   available germline subjects; [MutationEngine$set_germline_subject()]
//'   to set the active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject dataframe
//' head(m_engine$get_active_germline(), 5)
    .method("get_active_germline", &MutationEngine::get_active_germline)

//' @name MutationEngine$set_germline_subject
//' @title Setting the germline subject
//' @description This method sets the germline subject.
//' @details The subject must be one among those reported by
//'   [MutationEngine$get_germline_subjects()].
//' @return Set the germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to get the
//'   available germline subjects; [MutationEngine$get_active_germline()]
//'   to get the active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # set the active germline subject dataframe
//' m_engine$set_germline_subject("NA18941")
    .method("set_germline_subject", &MutationEngine::set_germline_subject)

//' @name MutationEngine$get_germline_subjects
//' @title Getting the germline subjects
//' @description This method returns the available germline subjects.
//' @details The germline subjects method returns a dataframe
//'   containing the available germline subjects. The column `sample`
//'   reports the subject name; the columns `pop` and `super_pop`
//'   contain the subject population and super population,
//'   respectively; the column `gender` declares the subject gender.
//' @return A dataframe the available germline subjects.
//' @seealso [MutationEngine$get_active_germline()] to get the
//'   available germline subjects; [MutationEngine$set_germline_subject()]
//'   to set the active germline.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject dataframe
//' head(m_engine$get_germline_subjects(), 5)
    .method("get_germline_subjects", &MutationEngine::get_germline_subjects)

//' @name MutationEngine$get_population_descriptions
//' @title Getting the population descriptions
//' @description This method returns the population descriptions.
//' @details The population descriptions are stored in a
//'   dataframe describing the populations. The column `code`
//'   contains the population codes; the columns `description`
//'   and `long description` report a brief and a long
//'   description for the populations, respectively.
//' @return A dataframe containing the population descriptions.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject dataframe
//' head(m_engine$get_population_descriptions(), 5)
    .method("get_population_descriptions",
            &MutationEngine::get_population_descriptions)

//' @name MutationEngine$get_species_rates
//' @title Getting the registered species rates
//' @description This method returns the registered species rates.
//' @details The registered species rates are returned in a
//'   dataframe. The column `species` contains the species names;
//'   the columns `SNV_rate`, `CNA_rate`, and `indel_rate`
//'   store the SNV, CNA, and indel rates, respectively.
//' @return A dataframe containing the registered species rates.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject dataframe
//' head(m_engine$get_species_rates(), 5)
    .method("get_species_rates",
            &MutationEngine::get_species_rates)

//' @name MutationEngine$get_SNV_signatures
//' @title Getting the SNV signatures
//' @description This method returns the available SNV
//'   signatures.
//' @details The signatures are returned in a dataframe
//'   containing the available SNV signatures and the
//'   corresponding mutation probability. The first column
//'   ("Type") describes a mutation in a context, while each
//'   of the remaining columns contains the probabilities
//'   of the mutations for one of the available SNV
//'   signatures.
//' @return A dataframe containing the available SNV
//'   signatures.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the indel dataframe
//' head(m_engine$get_SNV_signatures(), 5)
    .method("get_SNV_signatures",
            &MutationEngine::get_SNV_signatures_dataframe,
            "Get the SNV signatures dataframe")

//' @name MutationEngine$get_indel_signatures
//' @title Getting the indel signatures
//' @description This method returns the available indel
//'   signatures.
//' @details The signatures are returned in a dataframe
//'   containing the available indel signatures together with
//'   the corresponding mutation probability. The first column
//'   ("Type") describes a mutation in a context, while each
//'   of the remaining columns contains the probabilities of
//'   the mutations for one of the available indel signatures.
//' @return A dataframe containing the available  indel
//'   signatures.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the indel dataframe
//' head(m_engine$get_indel_signatures(), 5)
    .method("get_indel_signatures",
            &MutationEngine::get_indel_signatures_dataframe,
            "Get the indel signatures dataframe")

//' @name MutationEngine$get_known_drivers
//' @title Getting the known driver mutations
//' @description This method returns the known driver
//'   mutations.
//' @details The mutation are returned in a dataframe reporting
//'   the known driver mutations together with their types,
//'   associated tumours, affected genes, and code name. The
//'   first three columns ("`chr`", "`from`", and "`to`")
//'   report the mutation chromosome, the initial position
//'   and the final position, respectively. The next three
//'   columns ("`ref`", "`alt`", and "`mutation_type`")
//'   describe the reference sequence, the altered sequence,
//'   and the type of the mutation. The last three columns
//'   ("`tumour_type`", "`driver_gene`", and "`driver_code`")
//'   detail the tumour type associated to the mutation, the
//'   affected gene, and the driver code, which can be used
//'   to specify the mutation when adding mutants to the
//'   mutation engine.
//' @return A dataframe containing the known driver.
//' @seealso [MutationEngine$add_mutant()]
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the known driver dataframe
//' head(m_engine$get_known_drivers(), 5)
    .method("get_known_drivers",
            &MutationEngine::get_known_driver_mutations,
            "Get the known driver dataframe")

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
//'   be obtained by using the function [get_mutation_engine_codes()].
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
//' @seealso [get_mutation_engine_codes()] provides a list of the supported
//'   set-up codes.
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
//'   the index (optional: default value is 100).
//' @param COSMIC_account A named list containing "email" and "password" of
//'   a valid COSMIC account (required to download mutational signatures
//'   from COSMIC site).
//' @param max_index_size The maximum size of an admitted indel and, as a
//'   consequence, the maximum size of a motif stored in the repeated
//'   sequence index (optional: default value is 50).
//' @param max_repetition_storage The maximum number of repetitions per type
//'   stored in the repeated sequence index (optional: default value is
//'   500000).
//' @param tumour_type The type of tumour. This is currently used to select
//'   the admissible passenger CNAs. If any passenger CNA in the dataset is
//'   admissible, use the the empty string `""` (optional: default value is
//'   `""`).
//' @param tumour_study The nationality code of the tumour study. This is
//'   used to select the admissible passenger CNAs. If any tumor study in
//'   the dataset is admissible, use the the empty string `""` (optional:
//'   default value is `""`).
//' @param avoid_homozygous_losses An optional Boolean flag to avoid
//'   homozygous losses. When set to `TRUE`, passenger CNAs will be
//'   exclusively applied to regions covered by two alleles at least.
//'   (default: TRUE).
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso [MutationEngine$get_germline_subjects()] to get the available
//'   germline subjects; [MutationEngine$set_germline_subject()] to set the
//'   active germline subject; [MutationEngine$get_active_germline()] to get
//'   the active germline subject.
//' @return A mutation engine object.
//' @examples
//' # set the reference and SBS URLs
//' reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/release-111/",
//'                         "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
//'                         "dna.chromosome.22.fa.gz")
//' sbs_url <- paste0("https://raw.githubusercontent.com/",
//'                   "caravagnalab/ProCESS/refs/heads/1.0/inst/extdata/",
//'                   "SBS_demo_signatures.txt")
//' indel_url <- paste0("https://raw.githubusercontent.com/",
//'                     "caravagnalab/ProCESS/refs/heads/1.0/inst/extdata/",
//'                     "indel_demo_signatures.txt")
//' drivers_url <- paste0("https://raw.githubusercontent.com/",
//'                       "caravagnalab/ProCESS/refs/heads/1.0/inst/extdata/",
//'                       "driver_mutations_hg19.csv")
//' passenger_CNAs_url <- paste0("https://raw.githubusercontent.com/",
//'                              "caravagnalab/ProCESS/refs/heads/1.0/inst/",
//'                              "extdata/passenger_CNAs_hg19.csv")
//' germline_url <- paste0("https://zenodo.org/records/13166780/files/",
//'                        "germline_data_demo.tar.gz")
//'
//' # build a mutation engine and save the required files into the
//' # directory "Test". The `drivers_url` parameter is optional, but
//' # it is suggested to avoid passenger mutations on driver loci.
//' m_engine <- MutationEngine(directory = "Test",
//'                            reference_src = reference_url,
//'                            SBS_signatures_src = sbs_url,
//'                            indel_signatures_src = indel_url,
//'                            drivers_src = drivers_url,
//'                            passenger_CNAs_src = passenger_CNAs_url,
//'                            germline_src = germline_url)
//'
//' # if the parameters of a mutation engine construction match those of a
//' # previous construction, then the corresponding reference sequence,
//' # the SBS file, and the previously built context index are loaded from
//' # the set-up directory avoiding further computations.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'                            drivers_url, passenger_CNAs_url, germline_url)
//'
//' # if the `context_sampling` parameter changes, a new context index is
//' # built, while neither the reference sequence nor the SBS file are
//' # downloaded again.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'                            drivers_url, passenger_CNAs_url, germline_url,
//'                            context_sampling = 50)
//'
//' # a further construction with the same parameters avoids both
//' # downloads and context index construction.
//' m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
//'                            drivers_url, passenger_CNAs_url, germline_url,
//'                            context_sampling = 50)
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
//'                            COSMIC_account = list(email = "foo@bar.org",
//'                                                  password = "********"))
//' m_engine
//'
//' # remove the "demo" directory
//' unlink("demo", recursive = TRUE)
//'
//' # In alternative, pre-download the mutational signatures and pass their
//' # paths to `MutationEngine()` as parameters.
//' m_engine <- MutationEngine(setup_code = "demo",
//'                            SBS_signatures_src = "Test/SBS_signatures.txt",
//'                            indel_signatures_src = "Test/indel_signatures.txt")
//' m_engine
//'
//' # remove the "Test" and "demo" directories
//' unlink("Test", recursive = TRUE)
//' unlink("demo", recursive = TRUE)
  function("MutationEngine", &MutationEngine::build_MutationEngine,
           List::create(_["directory"] = "",
                        _["reference_src"] = "", _["SBS_signatures_src"] = "",
                        _["indel_signatures_src"] = "",
                        _["drivers_src"] = "", _["passenger_CNAs_src"] = "",
                        _["germline_src"] = "", _["setup_code"] = "",
                        _["COSMIC_account"] = R_NilValue,
                        _["germline_subject"] = "", _["context_sampling"] = 100,
                        _["max_motif_size"] = 50,
                        _["max_repetition_storage"] = 500000,
                        _["tumour_type"] = "", _["tumour_study"] = "",
                        _["avoid_homozygous_losses"] = true,
                        _["quiet"] = false),
           "Create a MutationEngine");

//' @name get_available_tumours_in
//' @title Getting the tumour types and studies available for a set-up code.
//' @description This method returns the tumour types and studies available
//'   for a set-up code.
//' @param setup_code The set-up code whose available tumour types and studies
//'   are requested.
//' @return A dataframe reporting the types and the studies available for a
//'   set-up code.
//' @seealso [MutationEngine()] to build a mutation engine
//' @export
//' @examples
//' # get the types and studies available for the "demo" set-up code
//' get_available_tumours_in("demo")
  function("get_available_tumours_in",
           &MutationEngine::get_available_tumour_type,
           List::create(_["setup_code"]),
           "Get the set of tumour types and studies for a set-up code.");

//' @name get_mutation_engine_codes
//' @title Getting the supported set-up codes
//' @description This method returns the supported codes for predefined set-up.
//' @return A dataframe reporting the code and a description for each
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
//'   This class is analogous to the class [SamplesForest],
//'   but each node is labelled with the mutations occurring
//'   for the first time on the cell represented by the node
//'   itself. Moreover each leaf is also associated with the
//'   genome mutations occurring in the corresponding cell.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'   cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A dataframe representing, for each of the identified
//'   cells, the identified (column `cell_id`), whenever the
//'   node is not a root, the ancestor identifier (column
//'   `ancestor`), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column `sample`), the mutant
//'   (column `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' }
//' @field get_first_occurrences Gets the identifier of the cell in which a
//'   mutation occurs for the first time\itemize{
//' \item \emph{Parameter:} \code{mutation} - A mutation being a
//'   SNV, a indel, or a CNA.
//' \item \emph{Return:} The identifier of the cell in which a mutation
//'   occurs for the first time.
//' }
//' @field get_germline_mutations Gets the germinal SNVs and indels\itemize{
//' \item \emph{Return:} A dataframe reporting `chr` (i.e., the
//'   chromosome), `chr_pos`" (i.e., the position in the chromosome),
//'   `allele` (in which the SNV occurs), `ref`, `alt`, `type` (i.e., either
//'   `"SNV"` or `"indel"`) and `class` (i.e., `"germinal"`).
//' }
//' @field get_germline_subject Gets the germline subject name\itemize{
//' \item \emph{Return:} The name of the subject whose germline is used.
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A dataframe representing, for each node
//'   in the forest, the identified (column `id`),
//'   whenever the node is not a root, the ancestor
//'   identifier (column `ancestor`), whenever the node
//'   was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column `sample`), the mutant (column
//'   `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' }
//' @field get_sampled_cell_CNAs Gets the CNAs of the sampled cells \itemize{
//' \item \emph{Returns:} A dataframe reporting `cell_id`, `type` (`"A"` for
//'   amplifications and `"D"` for deletions), `chr`, `begin`
//'   (i.e., the first CNA locus in the chromosome), `end` (i.e., the
//'   last CNA locus in the chromosome), `allele`, `src allele`
//'   (the allele origin for amplifications, `NA` for deletions), and
//'   `class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
//'   `"preneoplastic"`).
//' }
//' @field get_sampled_cell_mutations Gets the SNVs and the indels of the
//'   sampled cells\itemize{
//' \item \emph{Returns:} A dataframe reporting `cell_id`, `chr`, (i.e., the
//'   mutation chromosome), `begin` (i.e., position in the chromosome),
//'   `allele` (in which the SNV occurs), `ref`, `alt`, `type` (i.e., either
//'   `"SNV"` or `"indel"`), `cause`, and `class` (i.e., `"driver"`,
//'   `"passenger"`, `"germinal"` or `"preneoplastic"`) for each mutation
//'   in the sampled cell genomes.
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A dataframe containing, for each sample collected
//'   during the simulation, the columns "`name`", "`time`", "`ymin`",
//'   "`xmin`", "`ymax`", "`xmax`", "`tumour_cells`", and
//'   "`tumour_cells_in_bbox`". The columns "`ymin`", "`xmin`", "`ymax`",
//'   "`xmax`" report the boundaries of the sample bounding box, while
//'   "`tumour_cells`" and "`tumour_cells_in_bbox`" are the number of
//'   tumour cells in the sample and in the bounding box, respectively.
//' }
//' @field get_driver_mutations Getting the driver mutations\itemize{
//' \item \emph{Returns:} A dataframe reporting `mutant`, `order`,
//'   `type`, `CNA_type`, `chr`, `start`, `end`, `ref`, `alt`,
//'   `allele`, and `src_allele` for each driver mutations.
//' }
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A dataframe reporting `mutant`, `epistate`,
//'   `SNV_rate`, `indel_rate`, and `CNA_rate` for each registered species.
//' }
//' @field get_sticks Compute the forest sticks \itemize{
//' \item \emph{Returns:} The list of the forest sticks. Each stick is represented as
//'   the list of cell identifiers labelling the nodes in the stick
//'   from the higher to the deeper in the forest.
//' }
//' @field get_subforest_for Build a sub-forest using as leaves some of the original samples \itemize{
//' \item \emph{Parameter:} \code{sample_names} - The names of the samples whose cells will be used
//'   as leaves of the new forest.
//' \item \emph{Returns:} A samples forest built on the samples mentioned in `sample_names`.
//' }
//' @field get_absolute_chromosome_positions Get the absolute chromosome positions \itemize{
//' \item \emph{Returns:} A dataframe reporting the name (column "`chr`"), the length
//'   (column "`length`"), the initial absolute position (column "`from`"),
//'   and the final absolute position (column "`to`") of each chromosome.
//' }
//' @field save Save a phylogenetic forest in a file \itemize{
//' \item \emph{Parameter:} \code{filename} - The path of the file in which the phylogenetic
//'   forest must be saved.
//' }
  class_<PhylogeneticForest>("PhylogeneticForest")

//' @name PhylogeneticForest$get_nodes
//' @title Getting the forest nodes
//' @description This method returns the nodes of the forest.
//' @return A dataframe representing, for each node
//'   in the forest, the identified (column `cell_id`),
//'   whenever the node is not a root, the ancestor
//'   identifier (column `ancestor`), whenever the
//'   node was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column `sample`), the mutant (column
//'   `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' @seealso [SamplesForest$get_nodes()] for usage examples
    .method("get_nodes", (List (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_nodes),
            "Get the nodes of the forest")

//' @name PhylogeneticForest$get_coalescent_cells
//' @title Retrieving most recent common ancestors
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
//' @return A dataframe representing, for each of the identified
//'   cells, the identified (column `cell_id`), whenever the
//'   node is not a root, the ancestor identifier (column
//'   `ancestor`), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column `sample`), the mutant
//'   (column `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' @seealso [SamplesForest$get_coalescent_cells()] for usage examples
    .method("get_coalescent_cells",
            (List (PhylogeneticForest::*)(const std::list<RACES::Mutants::CellId>&) const)
                (&PhylogeneticForest::get_coalescent_cells),
            "Get the most recent common ancestor of some cells")
    .method("get_coalescent_cells",
            (List (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_coalescent_cells),
            "Get the most recent common ancestor of all the forest trees")

//' @name PhylogeneticForest$get_subforest_for
//' @title Building sub-forests
//' @description This method builds a sub-forest using as leaves some
//'   of the original samples.
//' @param sample_names The names of the samples whose cells will be used
//'   as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @seealso [SamplesForest$get_subforest_for()] for usage examples
    .method("get_subforest_for", &PhylogeneticForest::get_subforest_for,
            "Get the sub-forest for some of the original samples")

//' @name PhylogeneticForest$get_samples_info
//' @title Retrieving sample information
//' @description This method retrieves information about
//'   the samples whose cells were used as leaves
//'   of the samples forest.
//' @return A dataframe containing, for each sample collected during the
//'   simulation, the columns "`name`", "`time`", "`id`", "`ymin`", "`xmin`",
//'    "`ymax`", "`ymax`", "`xmax`", "`tumour_cells`", "`tumour_cells_in_bbox`",
//'   "`DNA_quantity`", and "`equivalent_normal_cells`". The columns "`ymin`",
//'   "`xmin`", "`ymax`", and "`xmax`" report the boundaries of the sample
//'   bounding box, while "`tumour_cells`" and "`tumour_cells_in_bbox`" are the
//'   number of tumour cells in the sample and in the bounding box,
//'   respectively. Finally, "`DNA_quantity`" contains the overall number of
//'   tumoral bases, i.e., the sum of the lengths of all the alleles of all the
//'   sample tumoral cells, and "`equivalent_normal_cells`" contains the number
//'   of normal cells that contain as much DNA as the sample tumoral cells.
//'   The "`DNA_quantity`" is stored as a real number despite being a natural
//'   number as it usually exceeds the largest natural number that can be
//'   natively represented by R.
//' @seealso [SamplesForest$get_samples_info()] for usage examples,
//'   [Simulation$get_samples_info()]
    .method("get_samples_info", &PhylogeneticForest::get_samples_info,
            "Get some pieces of information about the samples")

//' @name PhylogeneticForest$get_driver_mutations
//' @title Getting the driver mutations
//' @description This method returns the applied driver mutations.
//' @return A dataframe consisting in eight columns `mutant`, `order`, `type`,
//'    `CNA_type`, `chr`, `start`, `end`, `ref`, `alt`, `allele`, and
//'    `allele_srd`. Each row in the dataframe reports one driver mutations.
//'    The fields `mutant` and `order` report the name of the mutant and the
//'    application order among the mutant driver mutations, respectively.
//'    The column `type` declares the mutation type and contains "`SID`",
//'    "`CNA`", or "`WGD`" when the mutation is an SNV/indel, a CNA, or
//'    a whole genome duplication, respectively. When the mutation is a CNA,
//'    the `CNA_type` can either be "`A`" (i.e., amplification) or "`D`"
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
    .method("get_driver_mutations", &PhylogeneticForest::get_driver_mutations,
            "Get the applied driver mutations")

//' @name PhylogeneticForest$get_species_info
//' @title Getting the species
//' @description This method describes the simulated species.
//' @return A dataframe reporting `mutant`, `epistate`, `SNV_rate`,
//'   `indel_rate`, and `CNA_rate` for each registered species.
    .method("get_species_info", &PhylogeneticForest::get_species_info,
            "Get the recorded species")

//' @name PhylogeneticForest$get_germline_subject
//' @title Getting the germline subject
//' @description This method returns a dataframe reporting the germline
//'   subject name (column "sample"), population (column "pop"),
//'   super-population (column "super_pop"), and gender (column "gender").
//' @return The name of the subject whose germline is used.
    .method("get_germline_subject",
            &PhylogeneticForest::get_germline_subject_df,
            "Get the germline subject")

//' @name PhylogeneticForest$get_sampled_cell_CNAs
//' @title Getting the sampled cells CNAs
//' @description This method returns the CNAs of the sample cells.
//' @details This method builds a dataframe representing all the CNAs
//'   in the cells sampled during the simulation and represented by
//'   the leaves of the phylogenetic forest.
//' @param cell_id The identifier of the cell whose CNAs are aimed (optional).
//' @return A dataframe reporting `cell_id`, `type` (`"A"` for amplifications
//'   and `"D"` for deletions), `chr`, `begin` (i.e., the first CNA
//'   locus in the chromosome), `end` (i.e., last CNA locus in the chromosome),
//'   `allele`, `src allele` (the allele origin for amplifications, `NA` for
//'   deletions), and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"`
//'   or `"preneoplastic"`).
//' @examples
//' sim <- SpatialSimulation()
//'
//' sim$add_mutant(name = "Clone_1",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.05, "-" = 0.01))
//' sim$place_cell("Clone_1+", 500, 500)
//' sim$run_up_to_size("Clone_1-", 1000)
//'
//' sim$sample_cells("Sample_A", c(475, 475), c(525, 525))
//' sample_forest <- sim$get_samples_forest()
//'
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' m_engine$add_mutant("Clone_1", list("+" = c(SNV = 1e-7, indel = 1e-8),
//'                                     "-" = c(SNV = 3e-7, CNA = 1e-11)),
//'                     list(SNV("22", 10510210, "C", allele = 1),
//'                          CNA("D", "22", 5010000, 200000,
//'                                    allele = 1)))
//'
//' m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))
//'
//' phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
//'
//' mutations <- phylo_forest$get_sampled_cell_CNAs()
//'
//' head(mutations)
//' @seealso `PhylogeneticForest$$get_sampled_cell_mutations()`
    .method("get_sampled_cell_CNAs", (List (PhylogeneticForest::*)(const RACES::Mutants::CellId&) const)
                (&PhylogeneticForest::get_sampled_cell_CNAs),
            "Get the CNAs of a sampled cell")
    .method("get_sampled_cell_CNAs", (List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_sampled_cell_CNAs),
            "Get the CNAs of all the sampled cells")

//' @name PhylogeneticForest$get_sampled_cell_mutations
//' @title Getting the sampled cells mutations
//' @description This method returns the mutations of the sample cells.
//' @details This method builds a dataframe representing all the SNV
//'   and the indel mutations in the cells sampled during the simulation
//'   and represented by the leaves of the phylogenetic forest.
//'   The dataframe also reports the allele in which the mutations occur to
//'   support double occurrences due to CNAs.
//' @param cell_id The identifier of the cell whose mutations are aimed
//'   (optional).
//' @return A dataframe reporting `cell_id`, `chr`, (i.e., the mutation
//'   chromosome), `chr_pos` (i.e., position in the chromosome), `allele`
//'   (in which the mutation occurs), `ref`, `alt`, `type` (i.e., either
//'   `"SNV"` or `"indel"`), `cause`, and `class` (i.e., `"driver"`,
//'   `"passenger"`, `"germinal"` or `"preneoplastic"`) for each mutation
//'   in the sampled cell genomes.
//' @examples
//' sim <- SpatialSimulation()
//'
//' sim$add_mutant(name = "Clone_1",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.05, "-" = 0.01))
//' sim$place_cell("Clone_1+", 500, 500)
//' sim$run_up_to_size("Clone_1-", 1000)
//'
//' sim$sample_cells("Sample_A", c(475, 475), c(525, 525))
//' sample_forest <- sim$get_samples_forest()
//'
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' m_engine$add_mutant("Clone_1", list("+" = c(SNV = 1e-7, indel = 1e-8),
//'                                     "-" = c(SNV = 3e-7, CNA = 1e-11)),
//'                     list(SNV("22", 10510210, "C", allele = 1),
//'                          CNA("D", "22", 5010000, 200000,
//'                                    allele = 1)))
//'
//' m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))
//'
//' phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
//'
//' mutations <- phylo_forest$get_sampled_cell_mutations()
//'
//' head(mutations)
//' @seealso `PhylogeneticForest$$get_sampled_cell_CNAs()`
    .method("get_sampled_cell_mutations",
            (List (PhylogeneticForest::*)(const RACES::Mutants::CellId&) const)
                (&PhylogeneticForest::get_sampled_cell_SIDs),
            "Get the SNVs and the indels of a sampled cell")
    .method("get_sampled_cell_mutations", (List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_sampled_cell_SIDs),
            "Get the SNVs and indels of all the sampled cells")

//' @name PhylogeneticForest$get_germline_mutations
//' @title Getting the germinal mutations
//' @description This method returns the forest SNVs and indels.
//' @details Its builds a dataframe representing all the germinal
//'   SNVs and indels of the cells represented in the phylogenetic forest.
//'   The dataframe also reports the allele in which the mutations occur to
//'   support double occurrences due to CNAs.
//' @return A dataframe reporting "`chr`", "`chr_pos`" (i.e., the position in
//'   the chromosome), "`allele`" (in which the mutation occurs), "`ref`",
//'   "`alt`", "`cause`", "`type`" (i.e., either `"SNV"` or `"indel"`) and
//'   "`class`" (i.e., `"germinal"`).
//' @seealso `vignette("mutations")` for usage examples
    .method("get_germline_mutations", &PhylogeneticForest::get_germline_SIDs,
            "Get the germinal SNVs and indels")

//' @name PhylogeneticForest$get_absolute_chromosome_positions
//' @title Getting the absolute chromosome positions
//' @description This method returns the absolute chromosome positions.
//' @details Its builds a dataframe reporting the name, the length, and the
//'   initial and final absolute positions of each chromosome in the
//'   reference genome.
//' @return A dataframe reporting the name (column "`chr`"), the length
//'   (column "`length`"), the initial absolute position (column "`from`"),
//'   and the final absolute position (column "`to`") of each chromosome.
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
//' @seealso [SamplesForest$get_sticks()] for usage examples
    .method("get_sticks", (std::list<std::list<RACES::Mutants::CellId>> (PhylogeneticForest::*)(const double) const)(&PhylogeneticForest::get_sticks),
            "Get the forest sticks")
    .method("get_sticks", (std::list<std::list<RACES::Mutants::CellId>> (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_sticks),
            "Get the forest sticks")

//' @name PhylogeneticForest$get_exposures
//' @title Getting the timed exposure dataframe
//' @description This method returns a dataframe representing the exposure
//'   evolution over time.
//' @return A dataframe reporting `time`, `signature`, `exposure` and,
//'   `type`.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_exposures", &PhylogeneticForest::get_timed_exposures,
            "Get the timed exposure dataframe")

//' @name PhylogeneticForest$get_bulk_allelic_fragmentation
//' @title Getting the bulk allelic fragmentation dataframe
//' @description This method returns a dataframe representing the bulk allelic
//'     fragmentation of the genome.
//' @param sample_name The name of the sample whose bulk allelic fragmentation
//'     is aimed (optional).
//' @return A dataframe reporting, for each genomic fragment and for all
//'     the allelic type on the genomic fragment, the chromosome (`chr`),
//'     the first base position (`begin`), the last base position (`end`),
//'     the number of copy of the major and minor alleles (`major` and
//'     `minor`, respectively), and the ratio between the number of cells
//'     exhibiting this allelic type and the total number of cells in the
//'     sample.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_bulk_allelic_fragmentation",
            (Rcpp::List (PhylogeneticForest::*)(const std::string&) const)
                (&PhylogeneticForest::get_bulk_allelic_fragmentation),
            "Get the bulk allelic fragmentation dataframe")
    .method("get_bulk_allelic_fragmentation",
            (Rcpp::List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_bulk_allelic_fragmentation),
            "Get the bulk allelic fragmentation dataframe")

//' @name PhylogeneticForest$get_cell_allelic_fragmentation
//' @title Getting the cell allelic fragmentation dataframe
//' @description This method returns a dataframe representing the allelic
//'     fragmentation of each sampled cell.
//' @return A dataframe reporting, for each cell, for each genomic fragment,
//'     and for all the allelic type on the genomic fragment, the cell
//'     identifier (`cell_id`), the chromosome (`chr`), the first base
//'     position (`begin`), the last base position (`end`), and the number
//'     of copy of the major and minor alleles (`major` and `minor`,
//'     respectively).
//' @seealso `vignette("mutations")` for usage examples
    .method("get_cell_allelic_fragmentation",
            (Rcpp::List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_cell_allelic_fragmentation),
            "Get the cell allelic fragmentation dataframe")

//' @name PhylogeneticForest$get_first_occurrences
//' @title Getting the mutation first occurrences
//' @description This method returns the identifier of the cell in which a mutation occurs
//'   for the first time.
//' @param mutation A mutation being a SNV, a indel, or a CNA.
//' @return The identifier of the cell in which a mutation occurs for the first time.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_first_occurrences", (Rcpp::List (PhylogeneticForest::*)(const SEXP&) const)
                (&PhylogeneticForest::get_first_occurrence),
            "Get the identifier of the cell in which the mutation occurs for the first time")

//' @name PhylogeneticForest$get_reference_path
//' @title Getting the reference genome path
//' @description This method returns the reference genome path.
//' @return The reference genome path.
//' @seealso [PhylogeneticForest$set_reference_path()]
    .method("get_reference_path", (std::string (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_reference_path),
            "Get the reference genome path")

//' @name PhylogeneticForest$set_reference_path
//' @title Setting the reference genome path
//' @description This method returns the reference genome path.
//' @return The reference genome path.
//' @seealso [PhylogeneticForest$get_reference_path()]
    .method("set_reference_path", (void (PhylogeneticForest::*)(const std::string))
                (&PhylogeneticForest::set_reference_path),
            "Set the reference genome path")

//' @name PhylogeneticForest$save
//' @title Saving a phylogenetic forest
//' @description This method saves a phylogenetic forest in a file.
//' @param filename The path of the file in which the phylogenetic.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//'   forest must be saved.
    .method("save",
            (void (PhylogeneticForest::*)(const std::string&, const bool) const)
            &PhylogeneticForest::save,
            "Save a phylogenetic forest")
    .method("save",
            (void (PhylogeneticForest::*)(const std::string&) const)
            &PhylogeneticForest::save,
            "Save a phylogenetic forest")

    // show
    .method("show", &PhylogeneticForest::show,
            "Describe the PhylogeneticForest");

//' @name load_phylogenetic_forest
//' @title Loading a phylogenetic forest
//' @description This method loads a phylogenetic forest from a file.
//' @param filename The path of the file from which the phylogenetic
//'   forest must be load.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @return The load phylogenetic forest
  function("load_phylogenetic_forest",
           (PhylogeneticForest (*)(const std::string&, const bool))
           &PhylogeneticForest::load,
           List::create(_["filename"] = "", _["quiet"]=false),
           "Recover a phylogenetic forest");
}
