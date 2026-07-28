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

//' @name Mutation_class
//' @title Either an SNV or an indel
//' @description A class to represents SNVs and indels
//' @details The objects of this class represent either an
//'   SNV or an indel and can be created by [SNV()] and [Mutation()].
//'   They provide the following methods:
//'   - <code>[get_alt()](Mutation-cash-get_alt.md)</code> returns
//'     the sequence as altered by the mutation.
//'   - <code>[get_cause()](Mutation-cash-get_cause.md)</code> returns
//'     the cause of the mutation.
//'   - <code>[get_chromosome()](Mutation-cash-get_chromosome.md)</code>
//'     returns the chromosome in which the mutation occurred.
//'   - <code>[get_dataframe()](Mutation-cash-get_dataframe.md)</code>
//'     returns a data frame representing the mutation.
//'   - <code>[get_position_in_chromosome()](Mutation-cash-get_position_in_chromosome.md)</code>
//'     returns the position in the chromosome of the mutation.
//'   - <code>[get_ref()](Mutation-cash-get_ref.md)</code> returns
//'     the sequence as it was before the mutation.
//' @keywords internal
//' @seealso [Mutation()]
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
//' @keywords internal
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso <code>[Mutation](Mutation_class.md)</code>
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
//' @seealso [Mutation()] for building SNVs and indels,
//'   <code>[Mutation](Mutation_class.md)</code>
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
//' @details This function generalizes the function [SNV()] by constructing
//'   SNVs and indels. However, it necessitates the specification of the
//'   reference sequence, whereas [SNV()] can infer it from the reference
//'   sequence itself.
//'
//'   Another distinction between this function and [SNV()] lies in the order
//'   of the `ref`-`alt` parameter: in [SNV()], the alt parameter precedes the
//'   optional ref parameter, while `Mutation()` adopts the reverse order.
//' @usage Mutation(chr, from, ref, alt, allele, cause)
//' @param chr The name of the chromosome in which the indel occurs.
//' @param from The position in the chromosome where the indel occurs.
//' @param ref The reference sequence.
//' @param alt The mutation altered sequence.
//' @param allele The allele in which the mutation must occur (optional).
//' @param cause The cause of the mutation (optional).
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
//' @seealso [SNV()] for SNV creation, <code>[Mutation](Mutation_class.md)</code>
    function("Mutation", &SIDMut::build_SID,
             List::create(_["chr"], _["from"], _["ref"], _["alt"],
                          _["allele"] = R_NilValue, _["cause"] = ""),
             "Create either an SNV or a indel");

//' @name CNA
//' @title Creating a CNA
//' @description This function creates a <code>[CNA](CNA_class.md)</code> object.
//' @usage CNA(type, chr, chr_pos, len, allele = NULL, src_allele = NULL)
//' @param type The CNA type: either `"A"` or `"D"` for amplification and
//'   deletion, respectively.
//' @param chr The name of the chromosome in which the CNA occurs.
//' @param from The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the CNA occurs. (optional)
//' @param src_allele The allele from which the region is amplified. (optional,
//'   for amplification only)
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
//' @seealso [Amplification()],
//'   [Deletion()], <code>[CNA](CNA_class.md)</code>
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
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100)
//'
//' cna
//' @seealso  [Deletion()], [CNA()], <code>[CNA](CNA_class.md)</code>
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
//' @examples
//' # create a deletion CNA
//' cna <- Deletion("Y", 40020, 200)
//'
//' cna
//' @seealso  [Amplification()], [CNA()], <code>[CNA](CNA_class.md)</code>
    function("Deletion", &CNA::build_deletion,
             List::create(_["chr"], _["from"], _["len"], _["allele"] = R_NilValue),
             "Create a CNA deletion");

//' @name CNA_class
//' @title A copy number alteration
//' @description A class to represent CNA
//' @details This class represents copy number alterations (CNA).
//'   The objects of this class are built by [CNA()], [Amplification()],
//'   and [Deletion()]. They provide the following methods:
//'   - <code>[get_allele()](CNA-cash-get_allele.md)</code> returns the
//'     identifier of the CNA allele.
//'   - <code>[get_dataframe()](CNA-cash-get_dataframe.md)</code> returns
//'     a data frame representing the CNA.
//'   - <code>[get_length()](CNA-cash-get_length.md)</code> returns the
//'     length of the CNA.
//'   - <code>[get_position_in_chromosome()](CNA-cash-get_length.md)</code>
//'     returns the position of the CNA in its chromosome.
//'   - <code>[get_src_allele()](CNA-cash-get_src_allele.md)</code> returns
//'     the identifier of the allele from which the new allele is copied when
//'     the CNA is an amplification.
//'
//' @keywords internal
//' @seealso [CNA()], [Amplification()], [Deletion()]
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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
//' @seealso <code>[CNA](CNA_class.md)</code>
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

//' @name MutationEngine_class
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
//'   The initialisation of a <code>[MutationEngine](MutationEngine_class.md)</code>
//'   object requires a reference sequence and the SBS and ID mutational
//'   signatures. An SBS index and a ID index of the reference sequence
//'   are then automatically built. This process may take time depending
//'   on the size of the reference sequence. Hence, the downloaded files
//'   together with the built indices are saved on the disk for subsequent
//'   <code>[MutationEngine](MutationEngine_class.md)</code> constructions.
//'
//'   The objects of this class are generated by the function
//'   [MutationEngine()] and provide the following methods and properties:
//'   - <code>[add_exposure()](MutationEngine-cash-add_exposure.md)</code>
//'     adds exposures to the mutation engine.
//'   - <code>[add_mutant()](MutationEngine-cash-add_mutant.md)</code> adds
//'     a mutant genomic specification to the mutation engine.
//'   - <code>[change_rates_from()](MutationEngine-cash-change_rates_from.md)</code>
//'     changes the passenger rates from a specified time.
//'   - <code>[get_SNV_signatures()](MutationEngine-cash-get_SNV_signatures.md)</code>
//'     returns the SNV signatures.
//'   - <code>[get_active_germline()](MutationEngine-cash-get_active_germline.md)</code>
//'     returns the active germline subject.
//'   - <code>[get_genome_info()](MutationEngine-cash-get_genome_info.md)</code>
//'     returns a data frame reporting reference genome information.
//'   - <code>[get_germline_subjects()](MutationEngine-cash-get_germline_subjects.md)</code>
//'     returns the available germline subjects.
//'   - <code>[get_indel_signatures()](MutationEngine-cash-get_indel_signatures.md)</code>
//'     returns the indel signatures.
//'   - <code>[get_known_drivers()](MutationEngine-cash-get_known_drivers.md)</code>
//'     returns the data frame of the known driver mutations.
//'   - <code>[get_population_descriptions()](MutationEngine-cash-get_population_descriptions.md)</code>
//'     returns the available population descriptions.
//'   - <code>[get_species_info()](MutationEngine-cash-get_species_info.md)</code>
//'     returns the registered species and their rates.
//'   - <code>[infinite_sites_model](MutationEngine-cash-infinite_sites_model.md)</code>
//'     is a Boolean property to enable/disable the infinite sites model.
//'   - <code>[place_mutations()](MutationEngine-cash-place_mutations.md)</code>
//'     places mutations on a <code>[SampleForest]</code> object and generates
//'     a <code>[PhylogeneticForest]</code> object.
//'   - <code>[set_germline_subject()](MutationEngine-cash-set_germline_subject.md)</code>
//'     sets the active germline subject.
//'
//' @keywords internal
//' @seealso [MutationEngine()]
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso [MutationEngine$change_rates_from()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @title Changing the passenger rates from a specified time
//' @description This method changes the passenger rates from a specified time.
//' @details This method changes the passenger rates from a specified time. The
//'   rates before the specified time and those of the unspecified epigenetic
//'   states are not affected.
//' @param time The time from which the passenger rates are set.
//' @param mutant_name The mutant name.
//' @param passenger_rates The list of the passenger rates whose names are the
//'   epigenetic states of the species or a single rate, if the mutant
//'   does not have an epigenetic state.
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
//' @seealso [MutationEngine$add_mutant()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_active_germline(), 5)
//' @seealso [MutationEngine$get_germline_subjects()],
//'   [MutationEngine$set_germline_subject()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
        .method("get_active_germline", &MutationEngine::get_active_germline)

//' @name MutationEngine$set_germline_subject
//' @title Setting the germline subject
//' @description This method sets the germline subject.
//' @details The subject must be one among those reported by
//'   [MutationEngine$get_germline_subjects()].
//' @return Set the germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # set the active germline subject data frame
//' m_engine$set_germline_subject("NA18941")
//' @seealso [MutationEngine$get_germline_subjects()],
//'   [MutationEngine$get_active_germline()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_germline_subjects(), 5)
//' @seealso [MutationEngine$get_active_germline()],
//'   [MutationEngine$set_germline_subject()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
        .method("get_population_descriptions",
                &MutationEngine::get_population_descriptions)

//' @name MutationEngine$get_species_info
//' @title Getting the registered species and their rates
//' @description This method returns the registered species and
//'   their rates.
//' @details The registered species and their rates during the
//'   simulation are returned in a data frame. The column
//'   `mutant` contains the mutant names; the columns `time`,
//'   `SNV_rate`, `indel_rate`, and `CNA_rate` store the time
//'   from which rates hold, and the corresponding the SNV,
//'   indel, and CNA rates, respectively. If the simulation has
//'   epigenetic states, then the data frame also contains the
//'   column `epistate` to store the species epigenetic state
//'   names.
//' @return A data frame containing the registered species rates.
//' @examples
//' # build a mutation engine
//' m_engine <- MutationEngine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_species_info(), 5)
//' @seealso [SampleForest$get_species_info()],
//'   [PhylogeneticForest$get_species_info()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso [SampleForest$get_indel_signatures()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso [SampleForest$get_SNV_signatures()],
//'   <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
        .method("get_known_drivers", &MutationEngine::get_known_driver_mutations,
                "Get the known driver data frame")

        .method("show", &MutationEngine::show);

//' @name MutationEngine
//' @title Creating a mutation engine
//' @description This function downloads and sets up the data
//'   requires by a mutation engine. Finally, it builds an object
//'   of the class <code>[MutationEngine](MutationEngine_class.md)</code>
//'   itself.
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
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
//' @export
//' @examples
//' # get the types available for the "demo" set-up code
//' get_available_tumours_in("demo")
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
//'   [MutationEngine()]
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
//' @seealso <code>[MutationEngine](MutationEngine_class.md)</code>
//'   [MutationEngine()]
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
//'
//'   The objects of this class are built by using
//'   [MutationEngine] objects (see `[vignette("mutations")]`)
//'   and provide the following methods and properties:
//'   - <code>[get_absolute_chromosome_positions()](PhylogeneticForest-cash-get_absolute_chromosome_positions.md)</code>
//'     returns the absolute chromosome positions.
//'   - <code>[get_bulk_allelic_fragmentation()](PhylogeneticForest-cash-get_bulk_allelic_fragmentation.md)</code>
//'     returns the genome bulk allelic fragmentation.
//'   - <code>[get_cell_allelic_fragmentation()](PhylogeneticForest-cash-get_cell_allelic_fragmentation.md)</code>
//'     returns the cell allelic fragmentation.
//'   - <code>[get_coalescent_cells()](PhylogeneticForest-cash-get_coalescent_cells.md)</code>
//'     returns the most recent common ancestors of the sampled cells.
//'   - <code>[get_driver_mutations()](PhylogeneticForest-cash-get_driver_mutations.md)</code>
//'     returns the data frame of the driver mutations.
//'   - <code>[get_exposures()](PhylogeneticForest-cash-get_exposures.md)</code>
//'     returns the timed exposure data frame.
//'   - <code>[get_first_occurrences()](PhylogeneticForest-cash-get_first_occurrences.md)</code>
//'     returns the identifier of the first cell in which a mutation emerged.
//'   - <code>[get_germline_mutations()](PhylogeneticForest-cash-get_germline_mutations.md)</code>
//'     returns the data frame of the forest germinal mutations.
//'   - <code>[get_germline_subject()](PhylogeneticForest-cash-get_germline_subject.md)</code>
//'     returns the germline subject of the phylogenetic forest.
//'   - <code>[get_mutation_statistics()](PhylogeneticForest-cash-get_mutation_statistics.md)</code>
//'     returns the statistics about mutations on each node.
//'   - <code>[get_node()](PhylogeneticForest-cash-get_node.md)</code>
//'     returns an object of type <code>[PhylogeneticForestNode]</code>.
//'   - <code>[get_nodes()](PhylogeneticForest-cash-get_nodes.md)</code>
//'     returns information about the nodes in the forest.
//'   - <code>[get_reference_path()](PhylogeneticForest-cash-get_reference_path.md)</code>
//'     returns the path to the reference genome FASTA file.
//'   - <code>[get_sampled_cell_CNAs()](PhylogeneticForest-cash-get_sampled_cell_CNAs.md)</code>
//'     returns the data frame of the CNAs in the sampled cells.
//'   - <code>[get_sampled_cell_mutations()](PhylogeneticForest-cash-get_sampled_cell_mutations.md)</code>
//'     returns the data frame of the SNVs and indels in the sampled cells.
//'   - <code>[get_samples_info()](PhylogeneticForest-cash-get_samples_info.md)</code>
//'     returns information about the samples generating the forest.
//'   - <code>[get_species_info()](PhylogeneticForest-cash-get_species_info.md)</code>
//'     returns information about the simulated species.
//'   - <code>[get_sticks()](PhylogeneticForest-cash-get_sticks.md)</code>
//'     computes the forest sticks.
//'   - <code>[get_subforest_for()](PhylogeneticForest-cash-get_subforest_for.md)</code>
//'     extracts a sub-forest for some of the samples.
//'   - <code>[partition_samples()](PhylogeneticForest-cash-partition_samples.md)</code>
//'     partitions the samples according to a labelling.
//'   - <code>[represents_cell()](PhylogeneticForest-cash-represents_cell.md)</code>
//'     tests whether a cell having a given identifier is represented by the
//'     forest.
//'   - <code>[save()](PhylogeneticForest-cash-save.md)</code>
//'     saves the forest.
//'   - <code>[set_reference_path()](PhylogeneticForest-cash-set_reference_path.md)</code>
//'     sets the path to the reference genome FASTA file.
//' @seealso <code>[SampleForest]</code>, [MutationEngine$place_mutations()]
    class_<PhylogeneticForest>("PhylogeneticForest")
        REGISTER_FOREST_COMMON_FIELD(PhylogeneticForest)
        .method("partition_samples",
                &PhylogeneticForest::partition_samples,
                "Partition each sample in the forest according to a labelling function")
        .method("get_driver_mutations", &PhylogeneticForest::get_driver_mutations,
                "Get the applied driver mutations")
        .method("get_germline_subject", &PhylogeneticForest::get_germline_subject_df,
                "Get the germline subject")
        .method("get_sampled_cell_CNAs", &PhylogeneticForest::get_sampled_cell_CNAs,
                "Get the CNAs of all sampled cells")
        .method("get_sampled_cell_mutations",
                (Rcpp::DataFrame (PhylogeneticForest::*)() const)(
                    &PhylogeneticForest::get_sampled_cell_SIDs),
                "Get the SNVs and indels of all the sampled cells")
        .method("get_sampled_cell_mutations",
                (Rcpp::DataFrame (PhylogeneticForest::*)(const bool) const)(
                    &PhylogeneticForest::get_sampled_cell_SIDs),
                "Get the SNVs and indels of all the sampled cells")
        .method("get_germline_mutations", &PhylogeneticForest::get_germline_SIDs,
                "Get the germinal SNVs and indels")
        .method("get_absolute_chromosome_positions",
                &PhylogeneticForest::get_absolute_chromosome_positions,
                "Get the absolute chromosome positions")
        .method("get_exposures", &PhylogeneticForest::get_timed_exposures,
                "Get the timed exposure data frame")
        .method("get_bulk_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)(const std::string &)
                     const)(&PhylogeneticForest::get_bulk_allelic_fragmentation),
                "Get the bulk allelic fragmentation data frame")
        .method("get_bulk_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_bulk_allelic_fragmentation),
                "Get the bulk allelic fragmentation data frame")
        .method("get_cell_allelic_fragmentation",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_cell_allelic_fragmentation),
                "Get the cell allelic fragmentation data frame")
        .method("get_first_occurrences",
                (Rcpp::List (PhylogeneticForest::*)(const SEXP &)
                     const)(&PhylogeneticForest::get_first_occurrence),
                "Get the identifier of the cell in which the mutation occurs for the "
                "first time")
        .method("get_reference_path",
                &PhylogeneticForest::get_reference_path,
                "Get the reference genome path")
        .method("get_mutation_statistics",
                (Rcpp::List (PhylogeneticForest::*)()
                     const)(&PhylogeneticForest::get_mutation_statistics),
                "Get the statistics about node mutations")
        .method("set_reference_path",
                &PhylogeneticForest::set_reference_path,
                "Set the reference genome path");

//' @name PhylogeneticForestNode
//' @title The node of a phylogenetic forest
//' @description This class represents the nodes of a phylogenetic forest.
//' @details This class represents the nodes of a phylogenetic forest. It does not
//'   have a user constructor. Its objects are produced by [get_node_tour()] and
//'   <code>[PhylogeneticForest$get_node()]</code>.
//'
//'   The objects of this class provide the following methods and properties:
//'   - <code>[cell_id](PhylogeneticForest-cash-cell_id.md)</code>
//'     represents the identifier of the associated cell.
//'   - <code>[parent](PhylogeneticForest-cash-parent.md)</code>
//'     represents the node's parent.
//'   - <code>[children](PhylogeneticForest-cash-children.md)</code>
//'     represents a list of the node's children.
//'   - <code>[is_root](PhylogeneticForest-cash-is_root.md)</code>
//'     is a Boolean flag that is `TRUE` if and only if the node is a forest
//'     root.
//'   - <code>[is_leaf](PhylogeneticForest-cash-is_leaf.md)</code>
//'     is a Boolean flag that is `TRUE` if and only if the node is a forest
//'     leaf.
//'   - <code>[sample_name](PhylogeneticForest-cash-sample_name.md)</code>
//'     is the name of the sample that collected the associated cell.
//'   - <code>[birth_time](PhylogeneticForest-cash-birth_time.md)</code>
//'     is the birth time of the cell associated to the node.
//'   - <code>[death_time](PhylogeneticForest-cash-death_time.md)</code>
//'     is the death time of the cell associated to the node.
//'   - <code>[life_span](PhylogeneticForest-cash-life_span.md)</code>
//'     is the life span of the cell associated to the node.
//'   - <code>[species_name](PhylogeneticForest-cash-species_name.md)</code>
//'     is the name of the associated cell's species.
//'   - <code>[epistate_name](PhylogeneticForest-cash-epistate_name.md)</code>
//'     is the name of the associated cell's epigenetic state.
//'   - <code>[mutant_name](PhylogeneticForest-cash-mutant_name.md)</code>
//'     is the name of the associated cell's mutant.
//'   - <code>[arising_mutations](PhylogeneticForest-cash-arising_mutations.md)</code>
//'     stores the mutations arising in the associated cell.
//'   - <code>[get_genome()](PhylogeneticForest-cash-get_genome.md)</code>
//'     returns the genome of the associated cell.
//' @seealso [get_node_tour()], <code>[PhylogeneticForestNodeTour]</code>,
//'   <code>[SampleForestNode]</code>,
//'   `vignette("node_tour")`
   class_<PhylogeneticForestNode>("PhylogeneticForestNode")
        REGISTER_NODE_COMMON_FIELDS(PhylogeneticForestNode)
        .property("arising_mutations", &PhylogeneticForestNode::arising_mutations,
            "The mutations occurring for the first time in the associated cell (Read-only)")
        .method("get_genome", &PhylogeneticForestNode::get_genome,
            "Get the genome of the associated cell");

//' @name PhylogeneticForestNodeTour
//' @title An iterator class over phylogenetic forest nodes
//' @description Iterators over phylogenetic forest nodes.
//' @details This class represents iterators over phylogenetic forest nodes.
//'   The objects of this class are built by [get_node_tour()] and provide
//'   the following methods and properties:
//'   - <code>[node](PhylogeneticForestNodeTour-cash-node.md)</code>
//'     is an object of the class <code>[SampleForestNode]</code> and
//'     represents the node pointed by the iterator.
//'   - <code>[label](PhylogeneticForestNodeTour-cash-label.md)</code>
//'     represents the label of the of the node pointed by the
//'     iterator. The presence of this field depends on the type of the
//'     [get_node_tour()]'s parameters used to create the tour object
//'     (OPTIONAL).
//'   - <code>[genome](PhylogeneticForestNodeTour-cash-genome.md)</code>
//'     is an object of the class <code>[GenomeMutations]</code> that
//'     represent the genome of the node pointed by the iterator
//'     (OPTIONAL).
//'   - <code>[step()](PhylogeneticForestNodeTour-cash-step.md)</code>
//'     moves the iterator to the next node in the tour.
//'   - <code>[done](PhylogeneticForestNodeTour-cash-done.md)</code>
//'     is a Boolean flag that is set to `TRUE` only when the tour ended.
//'
//' @seealso [get_node_tour()], <code>[PhylogeneticForestNode]</code>,
//'   <code>[PhylogeneticForestNodeTour]</code>,
//'   `vignette("node_tour")`
    class_<PhylogeneticForestNodeTour>("PhylogeneticForestNodeTour")
        REGISTER_TOUR_COMMON_FIELDS(PhylogeneticForestNodeTour)
        .property("genome", &PhylogeneticForestNodeTour::get_genome,
            "The genome of the node pointed by the iterator");

//' @name GenomeFragment
//' @title Representing a genome fragment
//' @description This class represents genome fragment.
//' @details This class represents genome fragment.
//'   The objects of this class are built by <code>[GenomeMutations]</code>.
//'   They provides the following properties and methods:
//'   - <code>[get_CIGAR()](GenomeFragment-cash-get_CIGAR.md)</code> returns
//'     the CIGAR code of the fragment with respect to the reference genome.
//'   - <code>[get_covered_reference_region()](GenomeFragment-cash-get_covered_reference_region.md)</code>
//'     returns the reference genome region covered by the fragment.
//'   - <code>[get_mutations()](GenomeFragment-cash-get_mutations.md)</code>
//'     returns the mutations occurring on the fragment.
//'   - <code>[sequence](GenomeFragment-cash-sequence.md)</code> stores
//'     the fragment DNA sequence.
//' @seealso <code>[GenomeMutations]</code>, [get_node_tour()],
//'   `vignette("node_tour")`
    class_<GenomeFragment>("GenomeFragment")

//' @name GenomeFragment$get_CIGAR
//' @title Getting the fragment CIGAR
//' @description This method returns the CIGAR code of the fragment with respect
//'   to the reference genome.
//' @return The CIGAR code of the fragment with respect to the reference genome.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get a fragment
//' fragment <- genome$get_fragment("22", 0, 16085625, 100)
//'
//' # get its CIGAR with respect to the reference genome
//' fragment$get_CIGAR()
//' @seealso <code>[GenomeFragment]</code>
        .method("get_CIGAR", &GenomeFragment::get_CIGAR,
                "Get the fragment CIGAR")

//' @name GenomeFragment$sequence
//' @title The fragment sequence
//' @description This property is the fragment sequence.
//' @return The fragment sequence.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get a fragment
//' fragment <- genome$get_fragment("22", 0, 16085625, 100)
//'
//' # get the fragment sequence
//' fragment$sequence
//' @seealso <code>[GenomeFragment]</code>
        .property("sequence", &GenomeFragment::get_sequence,
                  "The fragment sequence")

//' @name GenomeFragment$get_mutations
//' @title Getting the mutations laying on a fragment
//' @description This method returns the data frame of the mutations
//'   laying on the fragment.
//' @return A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`, `alt`,
//'   `cause`, and `nature`. Each row represent a SID.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get a fragment
//' fragment <- genome$get_fragment("22", 0, 16085625, 100)
//'
//' # get the fragment mutations
//' fragment$get_mutations()
//' @seealso <code>[GenomeFragment]</code>
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
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get a fragment
//' fragment <- genome$get_fragment("22", 0, 16085625, 100)
//'
//' # get the covered reference region
//' fragment$get_covered_reference_region()
//' @seealso <code>[GenomeFragment]</code>
        .method("get_covered_reference_region",
                &GenomeFragment::get_covered_reference_region,
                "Get the covered region")

        .method("show", &GenomeFragment::show);

//' @name GenomeMutations
//' @title Representing cell genome
//' @description This class represents genome mutations of phylogenetic
//'   forest's cells.
//' @details The objects of this class are built by [get_node_tour()] and
//'   [PhylogeneticForestNode$get_genome()]. The provide the following
//'   methods:
//'   - <code>[get_CNAs()](GenomeMutations-cash-get_CNAs.md)</code>
//'     returns a data frame of the genome CNAs.
//'   - <code>[get_allele_fragments()](GenomeMutations-cash-get_allele_fragments.md)</code>
//'     returns genome allele fragments.
//'   - <code>[get_alleles_covering_ref_region()](GenomeMutations-cash-get_alleles_covering_ref_region.md)</code>
//'     returns the alleles covering a reference region.
//'   - <code>[get_fragment()](GenomeMutations-cash-get_fragment.md)</code>
//'     returns a fragment of the genome.
//'   - <code>[get_mutations()](GenomeMutations-cash-get_mutations.md)</code>
//'     returns a data frame of the genome mutations.
//'   - <code>[get_region_aligned_to_ref()](GenomeMutations-cash-get_region_aligned_to_ref.md)</code>
//'     returns information about the fragment aligning to the reference.
//'
//' @seealso [get_node_tour()],
//'   <code>[PhylogeneticForestNodeTour]</code>,
//'   `vignette("node_tour")`
    class_<GenomeMutations>("GenomeMutations")

//' @name GenomeMutations$get_mutations
//' @title Getting genome mutations
//' @description This method returns a data frame representing the SIDs in the genome.
//' @param with_germline A Boolean flag to add germline mutations (default: `TRUE`).
//' @return A data frame consisting of 7 columns: `chr`, `allele`, `from`, `ref`, `alt`,
//'   `cause`, and `nature`. Each row represent a SID.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get the first mutations in the genome
//' head(genome$get_mutations())
//' @seealso <code>[GenomeMutations]</code>
        .method("get_mutations",
            (Rcpp::DataFrame (GenomeMutations::*)() const)(&GenomeMutations::get_mutations),
            "Get the genome mutations")
        .method("get_mutations",
            (Rcpp::DataFrame (GenomeMutations::*)(const bool) const)(&GenomeMutations::get_mutations),
            "Get the genome mutations")

//' @name GenomeMutations$get_CNAs
//' @title Getting genome CNAs
//' @description This method returns a data frame representing the CNAs in the genome.
//' @return A data frame consisting of 8 columns: `chr`, `begin`, `end`, `type`,
//'   `allele`, `src.allele`, `cause`, and `nature`. Each row represent a SID.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get the first CNAs in the genome
//' head(genome$get_CNAs())
//' @seealso <code>[GenomeMutations]</code>
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
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get the genome allele fragments
//' genome$get_allele_fragments()
//' @seealso <code>[GenomeMutations]</code>
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
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get the region in chromosome 22 allele 0 aligning to the
//' # reference region from position 16085625 whose size is 100
//' genome$get_region_aligned_to_ref("22", 0, 16085625, 100)
//' @seealso <code>[GenomeMutations]</code>
        .method("get_region_aligned_to_ref",
                &GenomeMutations::region_aligning_on_reference,
            "Get information about the fragment aligning to the specified "
            "reference region")

//' @name GenomeMutations$get_alleles_covering_ref_region
//' @title Getting the alleles covering a reference region
//' @description This method returns the identifiers of the alleles that
//'   containing the specified region of the reference genome
//' @param chromosome_name The name of the chromosome of the reference region.
//' @param from The first position in the chromosome of the reference region.
//' @param size The size of the reference region.
//' @return A list of allele identifiers. Each identifier in the list
//'   corresponds to an allele containing the specified reference region.
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # the genome has 6 alleles
//' genome
//'
//' # get the alleles in the genome which covers the specified region
//' genome$get_alleles_covering_ref_region("22", 16085625, 100)
//' @seealso <code>[GenomeMutations]</code>
        .method("get_alleles_covering_ref_region",
                &GenomeMutations::get_alleles_covering_ref_region,
            "Get the alleles covering a region of the reference genome")

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
//' @examples
//' # use a phylogenetic forest example
//' forest <- example("PhylogeneticForest")
//'
//' # get the genome of the cell having 2 as the identifier
//' genome <- forest$get_node(2)$get_genome()
//'
//' # get a fragment
//' genome$get_fragment("22", 0, 16085625, 100)
//' @seealso <code>[GenomeFragment]</code>
//'   <code>[GenomeMutations$get_region_aligned_to_ref()]</code>,
//'   <code>[GenomeMutations]</code>
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
//' @usage load_phylogenetic_forest(filename, quiet)
//' @param filename The path of the file from which the phylogenetic
//'   forest must be loaded.
//' @param quiet An optional Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @return The loaded phylogenetic forest
//' @seealso [PhylogeneticForest$save()], [load_sample_forest()]
//'   [load_forest()]
    function("load_phylogenetic_forest",
             (PhylogeneticForest (*)(const std::string &,
                                     const bool))&PhylogeneticForest::load,
             List::create(_["filename"] = "", _["quiet"] = false),
             "Load a phylogenetic forest");
}
