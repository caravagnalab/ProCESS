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

#include <utils.hpp>

#include "cna.hpp"
#include "phylogenetic_forest.hpp"
#include "sid.hpp"
#include "simulation.hpp"
#include "mutation_engine.hpp"

#include "utility.hpp"

PhylogeneticForest::PhylogeneticForest() : CLONES::Mutations::PhylogeneticForest() {}

PhylogeneticForest::PhylogeneticForest(
    const CLONES::Mutations::PhylogeneticForest &orig,
    const GermlineSubject &germline_subject, const std::filesystem::path reference_path,
    const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
    const TimedMutationalExposure &timed_SBS_exposures,
    const TimedMutationalExposure &timed_indel_exposures)
    : CLONES::Mutations::PhylogeneticForest{orig}, germline_subject{germline_subject},
      reference_path{reference_path}, driver_codes{driver_codes}
{
    using namespace CLONES::Mutations;
    timed_exposures[MutationType::Type::SBS] = timed_SBS_exposures;
    timed_exposures[MutationType::Type::INDEL] = timed_indel_exposures;
}

PhylogeneticForest::PhylogeneticForest(
    CLONES::Mutations::PhylogeneticForest &&orig, const GermlineSubject &germline_subject,
    const std::filesystem::path reference_path,
    const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
    const TimedMutationalExposure &timed_SBS_exposures,
    const TimedMutationalExposure &timed_indel_exposures)
    : CLONES::Mutations::PhylogeneticForest{std::move(orig)},
      germline_subject{germline_subject}, reference_path{reference_path},
      driver_codes{driver_codes}
{
    using namespace CLONES::Mutations;
    timed_exposures[MutationType::Type::SBS] = timed_SBS_exposures;
    timed_exposures[MutationType::Type::INDEL] = timed_indel_exposures;
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> roots;

    for (const auto node : CLONES::Mutations::PhylogeneticForest::get_roots()) {
        roots.emplace_back(this, node);
    }

    return roots;
}

std::map<CLONES::Mutants::CellId, size_t>
PhylogeneticForest::get_total_mutations(const std::map<CLONES::Mutants::CellId, size_t>& new_mutations) const
{
    using namespace CLONES::Mutants;
    std::map<CellId, size_t> total_mutations;

    const auto& root_ids = get_root_cell_ids();

    std::list<CellId> next{root_ids.begin(), root_ids.end()};

    for (const auto& root: root_ids) {
        total_mutations.emplace(root, new_mutations.at(root));

        next.push_back(root);
    }

    while (!next.empty()) {
        const auto next_node = get_node(next.front());
        const auto front_muts = total_mutations.at(next.front());
        for (const auto child: next_node.children()) {
            const auto child_id = child.get_id();
            next.push_back(child_id);

            total_mutations.emplace(child_id, front_muts+new_mutations.at(child_id));
        }

        next.pop_front();
    }

    return total_mutations;
}

Rcpp::IntegerVector fill_vector_by_value(const std::map<CLONES::Mutants::CellId, size_t>& values)
{
    Rcpp::IntegerVector vector(values.size());

    size_t i{0};
    for (const auto& [cell_id, value] : values) {
        vector[i] = value;

        i++;
    }

    return vector;
}

Rcpp::List PhylogeneticForest::get_mutation_statistics() const
{
    const auto new_SID_maps = get_new_mutations<CLONES::Mutations::SID>(get_mutation_first_cells());
    const auto total_SID_maps = get_total_mutations(new_SID_maps);

    const auto new_SIDs = fill_vector_by_value(new_SID_maps);
    const auto total_SIDs = fill_vector_by_value(total_SID_maps);

    const auto new_CNA_maps = get_new_mutations<CLONES::Mutations::CNA>(get_CNA_first_cells());
    const auto total_CNA_maps = get_total_mutations(new_CNA_maps);

    const auto new_CNAs = fill_vector_by_value(new_CNA_maps);
    const auto total_CNAs = fill_vector_by_value(total_CNA_maps);

    using namespace Rcpp;

    IntegerVector ids(new_SID_maps.size());
    size_t i{0};
    for (const auto& [cell_id, mut]: new_SID_maps) {
        ids[i] = cell_id;

        ++i;
    }

    return DataFrame::create(_["cell_id"] = ids, _["new_SIDs"] = new_SIDs,
                             _["new_CNAs"] = new_CNAs, _["total_SIDs"] = total_SIDs,
                             _["total_CNAs"] = total_CNAs);
}

Rcpp::List PhylogeneticForest::get_samples_info() const
{
    auto info = TissueSimulation::get_samples_info(get_samples());

    const auto &samples = get_samples();
    std::vector<size_t> DNA(samples.size(), 0);
    for (const auto& [leaf_id, leaf_mutations] : get_leaf_mutation_tour(*this)) {
        DNA[get_coming_from().at(leaf_id)] += leaf_mutations.allelic_size();
    }

    std::map<std::string, size_t> sample_name_map;
    for (size_t j = 0; j < samples.size(); ++j) {
        sample_name_map[samples[j].get_name()] = j;
    }

    using namespace Rcpp;

    NumericVector DNA_quantities(DNA.size());
    NumericVector equivalent_normal_cells(DNA.size());
    StringVector name_col = info["name"];
    for (auto i = 0; i < name_col.size(); ++i) {
        const std::string sample_name = Rcpp::as<std::string>(name_col[i]);

        const auto j = sample_name_map.at(sample_name);

        DNA_quantities[j] = DNA[i];
        equivalent_normal_cells[j] = static_cast<double>(DNA[i])
                                / get_germline_mutations().allelic_size();
    }

    return DataFrame::create(_["name"] = info["name"], _["id"] = info["id"],
                             _["xmin"] = info["xmin"], _["ymin"] = info["ymin"],
                             _["xmax"] = info["xmax"], _["ymax"] = info["ymax"],
                             _["tumour_cells"] = info["tumour_cells"],
                             _["tumour_cells_in_bbox"] = info["tumour_cells_in_bbox"],
                             _["time"] = info["time"], _["DNA_quantity"] = DNA_quantities,
                             _["equivalent_normal_cells"] = equivalent_normal_cells);
}

std::string find_code(const CLONES::Mutations::MutationSpec<CLONES::Mutations::SID> &sid,
                      const std::map<CLONES::Mutations::SID, std::string> &driver_codes)
{
    const auto found = driver_codes.find(static_cast<CLONES::Mutations::SID>(sid));
    if (found == driver_codes.end()) {
        return "";
    }

    return (found->second);
}

inline void fill_SID_row(const CLONES::Mutations::MutationSpec<CLONES::Mutations::SID> &sid,
                         const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
                         Rcpp::CharacterVector &types, Rcpp::CharacterVector &CNA_types,
                         Rcpp::CharacterVector &chrs, Rcpp::CharacterVector &refs,
                         Rcpp::CharacterVector &alts, Rcpp::IntegerVector &starts,
                         Rcpp::IntegerVector &ends, Rcpp::IntegerVector &alleles,
                         Rcpp::IntegerVector &src_alleles,
                         Rcpp::CharacterVector &natures, Rcpp::CharacterVector &causes,
                         Rcpp::CharacterVector &codes, const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    types[i] = "SID";
    CNA_types[i] = NA_STRING;
    refs[i] = sid.ref;
    alts[i] = sid.alt;
    chrs[i] = GenomicPosition::chrtos(sid.chr_id);
    starts[i] = sid.position;
    ends[i] = sid.position + sid.ref.size() - 1;
    alleles[i] = sid.allele_id;
    src_alleles[i] = NA_INTEGER;
    natures[i] = sid.get_nature_description();
    causes[i] = sid.cause;

    const auto code = find_code(sid, driver_codes);

    if (code == "") {
        codes[i] = NA_STRING;
    } else {
        codes[i] = code;
    }
}

inline void fill_CNA_row(const CLONES::Mutations::CNA &cna, Rcpp::CharacterVector &types,
                         Rcpp::CharacterVector &CNA_types, Rcpp::CharacterVector &chrs,
                         Rcpp::CharacterVector &refs, Rcpp::CharacterVector &alts,
                         Rcpp::IntegerVector &starts, Rcpp::IntegerVector &ends,
                         Rcpp::IntegerVector &alleles, Rcpp::IntegerVector &src_alleles,
                         Rcpp::CharacterVector &natures, Rcpp::CharacterVector &causes,
                         Rcpp::CharacterVector &codes, const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    types[i] = "CNA";
    refs[i] = NA_STRING;
    alts[i] = NA_STRING;
    chrs[i] = GenomicPosition::chrtos(cna.chr_id);
    starts[i] = cna.position;
    ends[i] = cna.position + cna.length - 1;
    alleles[i] = (cna.dest == RANDOM_ALLELE ? NA_INTEGER : cna.dest);
    natures[i] = cna.get_nature_description();
    causes[i] = cna.cause;

    if (cna.type == CLONES::Mutations::CNA::Type::AMPLIFICATION) {
        src_alleles[i] = cna.source;
        CNA_types[i] = "A";
    } else {
        src_alleles[i] = NA_INTEGER;
        CNA_types[i] = "D";
    }
    codes[i] = NA_STRING;
}

inline void fill_WGD_row(const std::string& mutant_name,
                         Rcpp::CharacterVector &types, Rcpp::CharacterVector &CNA_types,
                         Rcpp::CharacterVector &chrs, Rcpp::CharacterVector &refs,
                         Rcpp::CharacterVector &alts, Rcpp::IntegerVector &starts,
                         Rcpp::IntegerVector &ends, Rcpp::IntegerVector &alleles,
                         Rcpp::IntegerVector &src_alleles,
                         Rcpp::CharacterVector &natures, Rcpp::CharacterVector &causes,
                         Rcpp::CharacterVector &codes, const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    types[i] = "WGD";
    CNA_types[i] = NA_STRING;
    refs[i] = NA_STRING;
    alts[i] = NA_STRING;
    chrs[i] = NA_INTEGER;
    starts[i] = NA_INTEGER;
    ends[i] = NA_INTEGER;
    alleles[i] = NA_INTEGER;
    src_alleles[i] = NA_INTEGER;
    natures[i] = CLONES::Mutations::Mutation::get_nature_description(CLONES::Mutations::Mutation::Nature::DRIVER);
    causes[i] = mutant_name;
    codes[i] = NA_STRING;
}

void fill_mutation_list(const CLONES::Mutations::MutationList& mutations,
                        const std::string& mutant_name,
                        const std::map<CLONES::Mutations::SID, std::string>& driver_codes,
                        Rcpp::CharacterVector &mutant_names, 
                        Rcpp::CharacterVector &types, Rcpp::CharacterVector &CNA_types,
                        Rcpp::CharacterVector &chrs, Rcpp::CharacterVector &refs,
                        Rcpp::CharacterVector &alts, Rcpp::IntegerVector &application_order,
                        Rcpp::IntegerVector &starts, Rcpp::IntegerVector &ends,
                        Rcpp::IntegerVector &alleles, Rcpp::IntegerVector &src_alleles,
                        Rcpp::CharacterVector &natures, Rcpp::CharacterVector &causes,
                        Rcpp::CharacterVector &codes, size_t &i, size_t &mut_i)
{
    using namespace CLONES::Mutations;

    for (auto dm_it = mutations.begin(); dm_it != mutations.end(); ++dm_it, ++mut_i, ++i) {
        mutant_names[i] = mutant_name;
        application_order[i] = mut_i;
        switch (dm_it.get_type()) {
        case MutationList::MutationType::SID_TURN:
            fill_SID_row(dm_it.get_last_SID(), driver_codes, types, CNA_types, chrs,
                         refs, alts, starts, ends, alleles, src_alleles, natures, causes,
                         codes, i);
            break;
        case MutationList::MutationType::CNA_TURN:
            fill_CNA_row(dm_it.get_last_CNA(), types, CNA_types, chrs, refs, alts,
                         starts, ends, alleles, src_alleles, natures, causes, codes, i);
            break;
        case MutationList::MutationType::WGD_TURN:
            fill_WGD_row(mutant_name, types, CNA_types, chrs, refs, alts, starts, ends,
                         alleles, src_alleles, natures, causes, codes, i);
            break;
        default:
            Rcpp::stop("Unsupported mutation type " +
                       std::to_string(dm_it.get_type()));
        }
    }
}

Rcpp::List PhylogeneticForest::const_node::arising_mutations() const
{
    using namespace Rcpp;

    const auto& mutations = node.arising_mutations();

    const size_t num_of_rows{mutations.size()};

    CharacterVector mutant_names(num_of_rows), chrs(num_of_rows), types(num_of_rows),
        CNA_types(num_of_rows), refs(num_of_rows), alts(num_of_rows), codes(num_of_rows),
        natures(num_of_rows), causes(num_of_rows);
    IntegerVector starts(num_of_rows), ends(num_of_rows), application_order(num_of_rows),
        alleles(num_of_rows), src_alleles(num_of_rows);

    size_t i{0}, mut_i{0};
    fill_mutation_list(mutations, node.get_species_name(), forest->driver_codes, 
                       mutant_names, types, CNA_types, chrs, refs, alts,
                       application_order, starts, ends, alleles, src_alleles, natures,
                       causes, codes, i, mut_i);

    return DataFrame::create(_["order"] = application_order,
                             _["type"] = types, _["CNA_type"] = CNA_types,
                             _["chr"] = chrs, _["start"] = starts, _["end"] = ends,
                             _["ref"] = refs, _["alt"] = alts, _["allele"] = alleles,
                             _["src_allele"] = src_alleles,
                             _["nature"] = natures, _["cause"] = causes,
                             _["code"] = codes);
}

Rcpp::List PhylogeneticForest::get_driver_mutations() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    const auto &mutational_properties = get_mutational_properties();

    size_t num_of_rows{0};
    for (const auto &[mutant, driver_mutations] :
         mutational_properties.get_driver_mutations()) {
        num_of_rows += driver_mutations.size();
    }

    CharacterVector mutant_names(num_of_rows), chrs(num_of_rows), types(num_of_rows),
        CNA_types(num_of_rows), refs(num_of_rows), alts(num_of_rows), codes(num_of_rows),
        natures(num_of_rows), causes(num_of_rows);
    IntegerVector starts(num_of_rows), ends(num_of_rows), application_order(num_of_rows),
        alleles(num_of_rows), src_alleles(num_of_rows);

    size_t i{0};
    for (const auto &[mutant, driver_mutations] :
         mutational_properties.get_driver_mutations()) {

        size_t mut_i{0};
        fill_mutation_list(driver_mutations, mutant, driver_codes,
                           mutant_names, types, CNA_types, chrs,
                           refs, alts, application_order, starts, ends, alleles,
                           src_alleles, natures, causes, codes, i, mut_i);
    }

    return DataFrame::create(_["mutant"] = mutant_names, _["order"] = application_order,
                             _["type"] = types, _["CNA_type"] = CNA_types,
                             _["chr"] = chrs, _["start"] = starts, _["end"] = ends,
                             _["ref"] = refs, _["alt"] = alts, _["allele"] = alleles,
                             _["src_allele"] = src_alleles, _["code"] = codes);
}

Rcpp::List PhylogeneticForest::get_species_info() const
{
    return MutationEngine::get_species_info(get_mutational_properties());
}

PhylogeneticForest
PhylogeneticForest::get_subforest_for(const std::vector<std::string> &sample_names) const
{
    PhylogeneticForest forest;

    forest.reference_path = reference_path;
    forest.timed_exposures = timed_exposures;

    static_cast<CLONES::Mutations::PhylogeneticForest &>(forest) =
        CLONES::Mutations::PhylogeneticForest::get_subforest_for(sample_names);

    return forest;
}

size_t count_mutations(const CLONES::Mutations::GenomeMutations &mutations)
{
    size_t counter{0};
    for (const auto &[chr_id, chromosome] : mutations.get_chromosomes()) {
        for (const auto &[allele_id, allele] : chromosome.get_alleles()) {
            for (const auto &[fragment_pos, fragment] : allele.get_fragments()) {
                counter += fragment.get_mutations().size();
            }
        }
    }

    return counter;
}

size_t
count_mutations(const std::map<CLONES::Mutants::CellId,
                               std::shared_ptr<CLONES::Mutations::CellGenomeMutations>>
                    &genome_mutations)
{
    size_t counter{0};
    for (const auto &[cell_id, mutations_ptr] : genome_mutations) {
        counter += count_mutations(*mutations_ptr);
    }

    return counter;
}

size_t
count_mutations(const CLONES::Mutations::MutationTour<CLONES::Mutations::GenomeMutations>& mutation_tour)
{
    size_t counter{0};

    for (const auto& [cell_id, node_mutations] : mutation_tour) {
        counter += count_mutations(node_mutations);
    }

    return counter;
}


size_t count_CNAs(const CLONES::Mutations::GenomeMutations &mutations)
{
    size_t counter{0};
    for (const auto &[chr_id, chromosome] : mutations.get_chromosomes()) {
        counter += chromosome.get_CNAs().size();
    }

    return counter;
}

size_t count_CNAs(const std::map<CLONES::Mutants::CellId,
                                 std::shared_ptr<CLONES::Mutations::CellGenomeMutations>>
                      &genome_mutations)
{
    size_t counter{0};
    for (const auto &[cell_id, mutations_ptr] : genome_mutations) {
        counter += count_CNAs(*mutations_ptr);
    }

    return counter;
}

size_t count_CNAs(const CLONES::Mutations::MutationTour<CLONES::Mutations::GenomeMutations>& mutation_tour)
{
    size_t counter{0};

    for (const auto& [cell_id, node_mutations] : mutation_tour) {
        counter += count_CNAs(node_mutations);
    }

    return counter;
}

void fill_mutation_lists(const CLONES::Mutants::CellId &cell_id,
                         const CLONES::Mutations::GenomeMutations &cell_mutations,
                         Rcpp::IntegerVector &cell_ids, Rcpp::CharacterVector &chr_names,
                         Rcpp::IntegerVector &chr_pos, Rcpp::IntegerVector &alleles,
                         Rcpp::CharacterVector &ref, Rcpp::CharacterVector &alt,
                         Rcpp::CharacterVector &types, Rcpp::CharacterVector &causes,
                         Rcpp::CharacterVector &classes, size_t &index)
{
    using namespace CLONES::Mutations;

    for (const auto &[chr_id, chromosome] : cell_mutations.get_chromosomes()) {
        for (const auto &[allele_id, allele] : chromosome.get_alleles()) {
            for (const auto &[fragment_pos, fragment] : allele.get_fragments()) {
                for (const auto &[mutation_pos, mutation_ptr] :
                     fragment.get_mutations()) {
                    cell_ids[index] = cell_id;
                    chr_names[index] = GenomicPosition::chrtos(chr_id);
                    chr_pos[index] = mutation_ptr->position;
                    alleles[index] = allele_id;
                    ref[index] = mutation_ptr->ref;
                    alt[index] = mutation_ptr->alt;
                    types[index] = (mutation_ptr->is_SBS() ? "SNV" : "indel");
                    causes[index] = mutation_ptr->cause;
                    classes[index] = mutation_ptr->get_nature_description();

                    ++index;
                }
            }
        }
    }
}

void fill_CNA_lists(const CLONES::Mutants::CellId &cell_id,
                    const CLONES::Mutations::GenomeMutations &cell_mutations,
                    Rcpp::IntegerVector &cell_ids, Rcpp::CharacterVector &chr_names,
                    Rcpp::IntegerVector &CNA_begins, Rcpp::IntegerVector &CNA_ends,
                    Rcpp::IntegerVector &src_alleles, Rcpp::IntegerVector &dst_alleles,
                    Rcpp::CharacterVector &types, Rcpp::CharacterVector &classes,
                    size_t &index)
{
    for (const auto &[chr_id, chromosome] : cell_mutations.get_chromosomes()) {
        for (const auto &cna_ptr : chromosome.get_CNAs()) {
            cell_ids[index] = cell_id;
            chr_names[index] = CLONES::Mutations::GenomicPosition::chrtos(chr_id);
            CNA_begins[index] = cna_ptr->begin();
            CNA_ends[index] = cna_ptr->end();
            bool is_amp = cna_ptr->type == CLONES::Mutations::CNA::Type::AMPLIFICATION;
            src_alleles[index] = (is_amp ? cna_ptr->source : NA_INTEGER);
            dst_alleles[index] = cna_ptr->dest;
            classes[index] = cna_ptr->get_nature_description();

            types[index] = (is_amp ? "A" : "D");

            ++index;
        }
    }
}

Rcpp::List PhylogeneticForest::get_absolute_chromosome_positions() const
{
    const auto &germline = get_germline_mutations();

    auto abspos = germline.get_absolute_chromosome_positions();

    using namespace Rcpp;

    NumericVector lengths(abspos.size()), froms(abspos.size()), tos(abspos.size());
    CharacterVector names(abspos.size());

    size_t index{0};
    size_t pos{1};
    for (const auto &[chr_id, from] : abspos) {
        names[index] = CLONES::Mutations::GenomicPosition::chrtos(chr_id);
        auto chr_length = germline.get_chromosome(chr_id).size();
        lengths[index] = chr_length;
        froms[index] = pos;
        pos += chr_length;
        tos[index] = (pos - 1);
        ++index;
    }

    return DataFrame::create(_["chr"] = names, _["length"] = lengths, _["from"] = froms,
                             _["to"] = tos);
}

Rcpp::List PhylogeneticForest::get_germline_SIDs() const
{
    const auto &germline = get_germline_mutations();
    size_t num_of_mutations{0};
    for (const auto &[chr_id, chr] : germline.get_chromosomes()) {
        for (const auto &[allele_id, allele] : chr.get_alleles()) {
            for (const auto &[f_pos, fragment] : allele.get_fragments()) {
                num_of_mutations += fragment.get_mutations().size();
            }
        }
    }

    using namespace Rcpp;

    IntegerVector chr_pos(num_of_mutations), alleles(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
        alt(num_of_mutations), types(num_of_mutations), causes(num_of_mutations),
        classes(num_of_mutations);

    size_t index{0};
    for (const auto &[chr_id, chr] : germline.get_chromosomes()) {
        for (const auto &[allele_id, allele] : chr.get_alleles()) {
            for (const auto &[f_pos, fragment] : allele.get_fragments()) {
                for (const auto &[mutation_pos, mutation_ptr] :
                     fragment.get_mutations()) {
                    chr_names[index] = CLONES::Mutations::GenomicPosition::chrtos(chr_id);
                    chr_pos[index] = mutation_ptr->position;
                    alleles[index] = allele_id;
                    ref[index] = mutation_ptr->ref;
                    alt[index] = mutation_ptr->alt;
                    types[index] = (mutation_ptr->is_SBS() ? "SNV" : "indel");
                    causes[index] = mutation_ptr->cause;
                    classes[index] = mutation_ptr->get_nature_description();

                    ++index;
                }
            }
        }
    }

    return DataFrame::create(_["chr"] = chr_names, _["from"] = chr_pos,
                             _["allele"] = alleles, _["ref"] = ref, _["alt"] = alt,
                             _["type"] = types, _["cause"] = causes,
                             _["class"] = classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_SIDs() const
{
    //size_t num_of_mutations = count_mutations(get_leaves_mutations());
    const auto leaf_tour = get_leaf_mutation_tour(*this);

    const size_t num_of_mutations = count_mutations(leaf_tour);

    using namespace Rcpp;

    IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
        alleles(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
        alt(num_of_mutations), types(num_of_mutations), causes(num_of_mutations),
        classes(num_of_mutations);

    size_t index{0};

    for (const auto& [leaf_id, leaf_mutations] : leaf_tour) {
        fill_mutation_lists(leaf_id, leaf_mutations, cell_ids, chr_names, chr_pos,
                            alleles, ref, alt, types, causes, classes, index);
    }

    return DataFrame::create(_["cell_id"] = cell_ids, _["chr"] = chr_names,
                             _["from"] = chr_pos, _["allele"] = alleles,
                             _["ref"] = ref, _["alt"] = alt, _["type"] = types,
                             _["cause"] = causes, _["class"] = classes);
}

Rcpp::List
PhylogeneticForest::get_cell_SIDs(const CLONES::Mutants::CellId &cell_id) const
{
    try {
        auto cell_mutations = get_cell_mutations(cell_id);

        return get_SID_dataframe(cell_mutations);
    } catch (std::domain_error& ex) {
         Rcpp::stop(ex.what());
    }
}

Rcpp::List PhylogeneticForest::get_SID_dataframe(
    const CLONES::Mutations::CellGenomeMutations &cell_mutations)
{
    const size_t num_of_mutations = count_mutations(cell_mutations);

    using namespace Rcpp;

    IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
        alleles(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
        alt(num_of_mutations), types(num_of_mutations), causes(num_of_mutations),
        classes(num_of_mutations);

    size_t index{0};
    fill_mutation_lists(cell_mutations.get_id(), cell_mutations, cell_ids,
                        chr_names, chr_pos, alleles, ref, alt, types, causes,
                        classes, index);

    return DataFrame::create(_["cell_id"] = cell_ids, _["chr"] = chr_names,
                             _["from"] = chr_pos, _["allele"] = alleles,
                             _["ref"] = ref, _["alt"] = alt, _["type"] = types,
                             _["cause"] = causes, _["class"] = classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_CNAs() const
{
    const auto leaf_tour = get_leaf_mutation_tour(*this);
    const size_t num_of_mutations = count_CNAs(leaf_tour);

    using namespace Rcpp;

    IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
        CNA_ends(num_of_mutations), src_alleles(num_of_mutations),
        dst_alleles(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), types(num_of_mutations),
        classes(num_of_mutations);

    size_t index{0};
    for (const auto& [leaf_id, leaf_mutations] : leaf_tour) {
        fill_CNA_lists(leaf_id, leaf_mutations, cell_ids, chr_names, CNA_begins,
                       CNA_ends, src_alleles, dst_alleles, types, classes, index);
    }

    return DataFrame::create(_["cell_id"] = cell_ids, _["type"] = types,
                             _["chr"] = chr_names, _["begin"] = CNA_begins,
                             _["end"] = CNA_ends, _["allele"] = dst_alleles,
                             _["src allele"] = src_alleles, _["class"] = classes);
}

Rcpp::List
PhylogeneticForest::get_cell_CNAs(const CLONES::Mutants::CellId &cell_id) const
{
    CLONES::Mutations::CellGenomeMutations cell_mutations;
    try {
        cell_mutations = get_cell_mutations(cell_id);
    } catch (std::domain_error& ex) {
        Rcpp::stop(ex.what());
    }

    size_t num_of_mutations = count_CNAs(cell_mutations);

    using namespace Rcpp;

    IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
        CNA_ends(num_of_mutations), src_alleles(num_of_mutations),
        dst_alleles(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), types(num_of_mutations),
        classes(num_of_mutations);

    size_t index{0};
    fill_CNA_lists(cell_id, cell_mutations, cell_ids, chr_names, CNA_begins,
                   CNA_ends, src_alleles, dst_alleles, types, classes, index);

    return DataFrame::create(_["cell_id"] = cell_ids, _["type"] = types,
                             _["chr"] = chr_names, _["begin"] = CNA_begins,
                             _["end"] = CNA_ends, _["allele"] = dst_alleles,
                             _["src allele"] = src_alleles, _["class"] = classes);
}

template <typename MUTATION_TYPE, typename R_MUTATION>
Rcpp::List get_first_occurrence(
    const std::map<MUTATION_TYPE, std::set<CLONES::Mutants::CellId>> &mutation_first_cells,
    const CLONES::Mutations::GenomeMutations &germline, const R_MUTATION &mutation)
{
    auto first_cell_it = mutation_first_cells.find(mutation);

    if (first_cell_it == mutation_first_cells.end()) {
        std::ostringstream oss;

        if constexpr (std::is_base_of_v<MUTATION_TYPE, CLONES::Mutations::SID>) {
            if (germline.includes(mutation)) {
                Rcpp::List R_cell_ids(1);

                R_cell_ids[0] = NA_INTEGER;
                return R_cell_ids;
            } else {
                oss << "The mutation " << mutation
                    << " does not occurs in the sampled cells.";
            }
        } else {
            oss << "The mutation " << mutation
                << " does not occurs in the sampled cells.";
        }

        Rcpp::stop(oss.str());
    }

    const auto &cell_ids = first_cell_it->second;

    Rcpp::List R_cell_ids(cell_ids.size());

    auto ids_it = cell_ids.begin();
    for (size_t i = 0; i < cell_ids.size(); ++i, ++ids_it) {
        R_cell_ids[i] = *ids_it;
    }
    return R_cell_ids;
}

Rcpp::List PhylogeneticForest::get_first_occurrence(const SEXP &mutation) const
{
    if (TYPEOF(mutation) == S4SXP) {

        Rcpp::S4 s4obj(mutation);
        if (s4obj.is("Rcpp_Mutation")) {
            Rcpp::Environment env(s4obj);
            Rcpp::XPtr<SIDMut> snv_ptr(env.get(".pointer"));

            return ::get_first_occurrence(get_mutation_first_cells(),
                                          get_germline_mutations(), *snv_ptr);
        }

        if (s4obj.is("Rcpp_CNA")) {
            Rcpp::Environment env(s4obj);
            Rcpp::XPtr<CNA> cna_ptr(env.get(".pointer"));

            return ::get_first_occurrence(get_CNA_first_cells(), get_germline_mutations(),
                                          *cna_ptr);
        }
    }

    ::Rf_error("This methods exclusively accepts as the parameters "
               "SNV, Indel, or CNA objects.");
}

void fill_timed_exposures(
    const std::map<CLONES::Time, CLONES::Mutations::MutationalExposure>
        &mutation_timed_exposures,
    const std::string &mutation_type_name, size_t &index, Rcpp::NumericVector &times,
    Rcpp::NumericVector &probs, Rcpp::CharacterVector &sig_names,
    Rcpp::CharacterVector &types)
{
    for (const auto &[time, exposure] : mutation_timed_exposures) {
        for (const auto &[sign_name, prob] : exposure) {
            times[index] = time;
            probs[index] = prob;
            sig_names[index] = sign_name;
            types[index] = mutation_type_name;
            ++index;
        }
    }
}

Rcpp::List PhylogeneticForest::get_timed_exposures() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;
    using namespace CLONES::Mutations;

    size_t dataframe_size{0};
    for (const auto &[type, mutation_timed_exposures] : timed_exposures) {
        (void)type;
        for (const auto &[time, exposure] : mutation_timed_exposures) {
            for (const auto &[sign_name, prob] : exposure) {
                (void)sign_name;
                ++dataframe_size;
            }
        }
    }

    NumericVector times(dataframe_size), probs(dataframe_size);
    CharacterVector sig_names(dataframe_size), types(dataframe_size);

    size_t index{0};
    fill_timed_exposures(timed_exposures.at(MutationType::Type::SBS), "SNV", index, times,
                         probs, sig_names, types);
    fill_timed_exposures(timed_exposures.at(MutationType::Type::INDEL), "indel", index,
                         times, probs, sig_names, types);

    return DataFrame::create(_["time"] = times, _["signature"] = sig_names,
                             _["exposure"] = probs, _["type"] = types);
}

size_t count_rows_in_allelic_bulk_data(
    const std::map<CLONES::Mutations::AllelicType, size_t> &chr_allelic_count,
    const double &num_of_cells)
{
    using namespace CLONES::Mutations;

    double total_count{0.0};
    size_t num_of_rows{0};

    for (const auto &[atype, count] : chr_allelic_count) {
        total_count += count;
        ++num_of_rows;
    }

    if (num_of_cells >= total_count + 1) {
        ++num_of_rows;
    }

    return num_of_rows;
}

size_t count_rows_in_bulk_allelic_fragmentation(
    const CLONES::Mutations::PhylogeneticForest::AllelicCount &allelic_count,
    const double &num_of_cells)
{
    size_t num_of_rows{0};
    for (const auto &[chr_id, chr_allelic_count] : allelic_count) {
        for (const auto &[pos, ca_count] : chr_allelic_count) {
            num_of_rows += count_rows_in_allelic_bulk_data(ca_count, num_of_cells);
        }
    }

    return num_of_rows;
}

void fill_allelic_bulk_data(
    const std::map<CLONES::Mutations::AllelicType, size_t> &chr_allelic_count,
    const CLONES::Mutations::ChromosomeId &chr_id,
    const CLONES::Mutations::ChrPosition frag_begin,
    const CLONES::Mutations::ChrPosition frag_end, const double &num_of_cells,
    Rcpp::StringVector &chromosomes, Rcpp::IntegerVector &fragment_begins,
    Rcpp::IntegerVector &fragment_ends, Rcpp::IntegerVector &major_counts,
    Rcpp::IntegerVector &minor_counts, Rcpp::NumericVector &ratios, size_t &row_idx)
{
    using namespace CLONES::Mutations;

    double total_count{0.0};

    for (const auto &[atype, count] : chr_allelic_count) {
        chromosomes[row_idx] = GenomicPosition::chrtos(chr_id);
        fragment_begins[row_idx] = frag_begin;
        fragment_ends[row_idx] = frag_end;

        if (atype[0] < atype[1]) {
            major_counts[row_idx] = atype[1];
            minor_counts[row_idx] = atype[0];
        } else {
            major_counts[row_idx] = atype[0];
            minor_counts[row_idx] = atype[1];
        }

        ratios[row_idx] = count / num_of_cells;
        total_count += count;

        ++row_idx;
    }

    if (num_of_cells >= total_count + 1) {
        chromosomes[row_idx] = GenomicPosition::chrtos(chr_id);
        fragment_begins[row_idx] = frag_begin;
        fragment_ends[row_idx] = frag_end;

        major_counts[row_idx] = 0;
        minor_counts[row_idx] = 0;

        ratios[row_idx] = (num_of_cells - total_count) / num_of_cells;

        ++row_idx;
    }
}

const std::list<CLONES::Mutants::CellId> &
PhylogeneticForest::get_cell_ids_in(const std::string &sample_name) const
{
    for (const auto &sample : get_samples()) {
        if (sample.get_name() == sample_name) {
            return sample.get_cell_ids();
        }
    }

    Rcpp::stop("Unknown sample \"" + sample_name + "\".");
}

Rcpp::List get_bulk_allelic_fragmentation(
    const CLONES::Mutations::PhylogeneticForest::AllelicCount &allelic_count,
    const std::map<CLONES::Mutations::ChromosomeId, CLONES::Mutations::ChromosomeMutations>
        &chr_map,
    const double &num_of_cells)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    const size_t num_of_rows =
        count_rows_in_bulk_allelic_fragmentation(allelic_count, num_of_cells);

    StringVector chromosomes(num_of_rows);
    IntegerVector fragment_begins(num_of_rows), fragment_ends(num_of_rows);
    IntegerVector major_counts(num_of_rows), minor_counts(num_of_rows);
    NumericVector ratios(num_of_rows);

    size_t row_idx{0};
    for (const auto &[chr_id, chr_allelic_count] : allelic_count) {
        auto a_count_it = chr_allelic_count.begin();

        if (a_count_it != chr_allelic_count.end()) {
            auto next_a_count_it = a_count_it;
            while (++next_a_count_it != chr_allelic_count.end()) {
                fill_allelic_bulk_data(a_count_it->second, chr_id, a_count_it->first,
                                       next_a_count_it->first - 1, num_of_cells,
                                       chromosomes, fragment_begins, fragment_ends,
                                       major_counts, minor_counts, ratios, row_idx);

                a_count_it = next_a_count_it;
            }

            fill_allelic_bulk_data(a_count_it->second, chr_id, a_count_it->first,
                                   chr_map.at(chr_id).size(), num_of_cells, chromosomes,
                                   fragment_begins, fragment_ends, major_counts,
                                   minor_counts, ratios, row_idx);
        }
    }

    return DataFrame::create(_["chr"] = chromosomes, _["begin"] = fragment_begins,
                             _["end"] = fragment_ends, _["major"] = major_counts,
                             _["minor"] = minor_counts, _["ratio"] = ratios);
}

Rcpp::List PhylogeneticForest::get_bulk_allelic_fragmentation() const
{
    const double num_of_cells = num_of_leaves();
    const auto &chr_map = get_germline_mutations().get_chromosomes();
    const auto allelic_count = get_allelic_count(2);

    return ::get_bulk_allelic_fragmentation(allelic_count, chr_map, num_of_cells);
}

Rcpp::List
PhylogeneticForest::get_bulk_allelic_fragmentation(const std::string &sample_name) const
{
    const double num_of_cells = get_cell_ids_in(sample_name).size();
    const auto &chr_map = get_germline_mutations().get_chromosomes();
    const auto allelic_count = get_allelic_count(sample_name, 2);

    return ::get_bulk_allelic_fragmentation(allelic_count, chr_map, num_of_cells);
}

void fill_allelic_cell_data(const CLONES::Mutations::AllelicType &allelic_type,
                            const CLONES::Mutants::CellId &cell_id,
                            const CLONES::Mutations::ChromosomeId &chr_id,
                            const CLONES::Mutations::ChrPosition frag_begin,
                            const CLONES::Mutations::ChrPosition frag_end,
                            Rcpp::IntegerVector &ids, Rcpp::StringVector &chromosomes,
                            Rcpp::IntegerVector &fragment_begins,
                            Rcpp::IntegerVector &fragment_ends,
                            Rcpp::IntegerVector &major_counts,
                            Rcpp::IntegerVector &minor_counts, size_t &row_idx)
{
    using namespace CLONES::Mutations;

    ids[row_idx] = cell_id;
    chromosomes[row_idx] = GenomicPosition::chrtos(chr_id);
    fragment_begins[row_idx] = frag_begin;
    fragment_ends[row_idx] = frag_end;

    if (allelic_type[0] < allelic_type[1]) {
        major_counts[row_idx] = allelic_type[1];
        minor_counts[row_idx] = allelic_type[0];
    } else {
        major_counts[row_idx] = allelic_type[0];
        minor_counts[row_idx] = allelic_type[1];
    }

    ++row_idx;
}

size_t count_rows_in_cell_allelic_fragmentation(
    const std::map<CLONES::Mutants::CellId,
                   std::shared_ptr<CLONES::Mutations::CellGenomeMutations>>
        &leave_mutations)
{
    size_t num_of_rows{0};

    for (const auto &[cell_id, mutations] : leave_mutations) {
        const auto b_points = mutations->get_CNA_break_points();
        const auto allelic_map = mutations->get_allelic_map(b_points, 2);

        for (const auto &[chr_id, chr_allelic_map] : allelic_map) {
            num_of_rows += chr_allelic_map.size();
        }
    }

    return num_of_rows;
}

size_t count_rows_in_cell_allelic_fragmentation(
    const CLONES::Mutations::MutationTour<CLONES::Mutations::GenomeMutations>& node_tour)
{
    size_t num_of_rows{0};

    for (const auto& [cell_id, cell_mutations] : node_tour) {
        const auto b_points = cell_mutations.get_CNA_break_points();
        const auto allelic_map = cell_mutations.get_allelic_map(b_points, 2);

        for (const auto &[chr_id, chr_allelic_map] : allelic_map) {
            num_of_rows += chr_allelic_map.size();
        }
    }

    return num_of_rows;
}

Rcpp::List PhylogeneticForest::get_cell_allelic_fragmentation() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    const auto leaf_tour = get_leaf_mutation_tour(*this);
    const size_t num_of_rows =
    //    count_rows_in_cell_allelic_fragmentation(get_leaves_mutations());
        count_rows_in_cell_allelic_fragmentation(leaf_tour);

    IntegerVector ids(num_of_rows);
    StringVector chromosomes(num_of_rows);
    IntegerVector fragment_begins(num_of_rows), fragment_ends(num_of_rows);
    IntegerVector major_counts(num_of_rows), minor_counts(num_of_rows);

    const auto &chr_map = get_germline_mutations().get_chromosomes();

    size_t row_idx{0};

    for (const auto& [leaf_id, leaf_mutations]: leaf_tour) {
        const auto b_points = leaf_mutations.get_CNA_break_points();
        const auto allelic_map = leaf_mutations.get_allelic_map(b_points, 2);

        for (const auto &[chr_id, chr_allelic_map] : allelic_map) {
            auto a_map_it = chr_allelic_map.begin();
            auto next_a_map_it = a_map_it;
            while (++next_a_map_it != chr_allelic_map.end()) {
                fill_allelic_cell_data(a_map_it->second, leaf_id, chr_id, a_map_it->first,
                                       next_a_map_it->first - 1, ids, chromosomes,
                                       fragment_begins, fragment_ends, major_counts,
                                       minor_counts, row_idx);

                a_map_it = next_a_map_it;
            }

            fill_allelic_cell_data(a_map_it->second, leaf_id, chr_id, a_map_it->first,
                                   chr_map.at(chr_id).size(), ids, chromosomes,
                                   fragment_begins, fragment_ends, major_counts,
                                   minor_counts, row_idx);
        }
    }

    return DataFrame::create(_["cell_id"] = ids, _["chr"] = chromosomes,
                             _["begin"] = fragment_begins, _["end"] = fragment_ends,
                             _["major"] = major_counts, _["minor"] = minor_counts);
}

void PhylogeneticForest::set_reference_path(const std::filesystem::path reference_path)
{
    if (!std::filesystem::exists(reference_path)) {
        Rcpp::stop("The reference genome file \"" + to_string(reference_path) +
                   "\" does not exists.");
    }

    this->reference_path = std::filesystem::absolute(reference_path);
}

void PhylogeneticForest::save(const std::string &filename, const bool quiet) const
{
    CLONES::Archive::Binary::Out out_archive(filename);

    auto ref_str = to_string(reference_path);

    out_archive & germline_subject & ref_str & driver_codes & timed_exposures;

    CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    out_archive.save(static_cast<const CLONES::Mutations::PhylogeneticForest &>(*this),
                     progress_bar, "forest");
}

PhylogeneticForest PhylogeneticForest::load(const std::string &filename, const bool quiet)
{
    if (!std::filesystem::exists(filename)) {
        Rcpp::stop("The file \"" + filename + "\" does not exist.");
    }

    if (!std::filesystem::is_regular_file(filename)) {
        Rcpp::stop("The file \"" + filename + "\" is not a regular file.");
    }

    CLONES::Archive::Binary::In in_archive(filename);

    PhylogeneticForest forest;

    in_archive & forest.germline_subject;

    std::string reference_path;
    in_archive & reference_path;
    forest.reference_path = reference_path;

    in_archive & forest.driver_codes & forest.timed_exposures;

    try {
        CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

        in_archive.load(static_cast<CLONES::Mutations::PhylogeneticForest &>(forest),
                        progress_bar, "forest");
    } catch (CLONES::Archive::WrongFileFormatDescr &ex) {
        raise_error(ex, "phylogenetic forest");
    } catch (CLONES::Archive::WrongFileFormatVersion &ex) {
        raise_error(ex, "phylogenetic forest");
    }

    return forest;
}

void PhylogeneticForest::show() const
{
    using namespace Rcpp;

    size_t num_of_leaves{0};
    for (const auto &sample : get_samples()) {
        num_of_leaves += sample.get_cell_ids().size();
    }

    Rcout << "PhylogeneticForest" << std::endl
          << "  # of trees: " << get_roots().size() << std::endl
          << "  # of nodes: " << num_of_nodes() << std::endl
          << "  # of leaves: " << num_of_leaves << std::endl
          << "  samples: {";

    std::string sep = "";
    for (const auto &sample : get_samples()) {
        Rcout << sep << "\"" << sample.get_name() << "\"";
        sep = ", ";
    }

    Rcout << "}" << std::endl
          << std::endl
          << "  # of emerged SNVs and indels: " << get_mutation_first_cells().size()
          << std::endl
          << "  # of emerged CNAs: " << get_CNA_first_cells().size() << std::endl
          << std::endl;
}
