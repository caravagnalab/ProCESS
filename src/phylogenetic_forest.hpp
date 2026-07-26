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

#ifndef __PROCESS_PHYLOGENETIC_FOREST__
#define __PROCESS_PHYLOGENETIC_FOREST__

#include <map>
#include <vector>

#include <Rcpp.h>

#include <mutation_engine.hpp>
#include <phylogenetic_forest.hpp>

#include "forest.hpp"
#include "genomic_data_storage.hpp"

class MutationEngine;

class PhylogeneticForest : public CLONES::Mutations::PhylogeneticForest
{
    GermlineSubject germline_subject;
    std::filesystem::path reference_path;

    std::map<CLONES::Mutations::SID, std::string> driver_codes;

    using TimedMutationalExposure =
        std::map<CLONES::Time, CLONES::Mutations::MutationalExposure>;

    std::map<CLONES::Mutations::MutationType::Type, TimedMutationalExposure>
        timed_exposures;

    PhylogeneticForest(const CLONES::Mutations::PhylogeneticForest &orig,
                       const GermlineSubject &germline_subject,
                       const std::filesystem::path reference_path,
                       const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
                       const TimedMutationalExposure &timed_SBS_exposures,
                       const TimedMutationalExposure &timed_indel_exposures);

    PhylogeneticForest(CLONES::Mutations::PhylogeneticForest &&orig,
                       const GermlineSubject &germline_subject,
                       const std::filesystem::path reference_path,
                       const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
                       const TimedMutationalExposure &timed_SBS_exposures,
                       const TimedMutationalExposure &timed_indel_exposures);

    template<typename MUTATION_TYPE>
    std::map<CLONES::Mutants::CellId, size_t>
    get_new_mutations(const std::map<MUTATION_TYPE, std::set<CLONES::Mutants::CellId>>& first_occurrences) const
    {
        std::map<CLONES::Mutants::CellId, size_t> new_mutations;

        for (const auto& [cell_id, cell] : get_cells()) {
            new_mutations.emplace(cell_id, 0);
        }

        for (const auto& [mutation, cell_ids]: first_occurrences) {
            for (const auto cell_id: cell_ids) {
                ++(new_mutations[cell_id]);
            }
        }

        return new_mutations;
    }

    std::map<CLONES::Mutants::CellId, size_t>
    get_total_mutations(const std::map<CLONES::Mutants::CellId, size_t>& new_mutations) const;

  public:
    inline static const std::string file_format_header = "ProCESS Phylogenetic Forest";
    inline static const uint8_t file_format_number = 0;

    using base_type = CLONES::Mutations::PhylogeneticForest;

    using const_node = ForestCore::const_node<PhylogeneticForest>;

    PhylogeneticForest();

    inline Rcpp::LogicalVector represents_cell(const SEXP cell_ids) const
    {
        return ForestCore::represents_cell(*this, cell_ids);
    }

    std::vector<const_node> get_roots() const;

    inline Rcpp::List get_nodes() const
    {
        return ForestCore::get_nodes(
            static_cast<const CLONES::Mutations::PhylogeneticForest &>(*this));
    }

    inline const_node get_node(const CLONES::Mutants::CellId &cell_id) const
    {
        return PhylogeneticForest::const_node(*this, cell_id);
    }

    inline const_node get_node(const SEXP cell_id) const
    {
        return PhylogeneticForest::const_node(*this, cell_id);
    }

    Rcpp::DataFrame get_mutation_statistics() const;

    const std::list<CLONES::Mutants::CellId> &
    get_cell_ids_in(const std::string &sample_name) const;

    Rcpp::DataFrame get_samples_info() const;

    Rcpp::DataFrame get_driver_mutations() const;

    Rcpp::List get_species_info() const;

    inline Rcpp::List get_coalescent_cells() const
    {
        return ForestCore::get_coalescent_cells(
            static_cast<const CLONES::Mutations::PhylogeneticForest &>(*this));
    }

    inline Rcpp::List
    get_coalescent_cells(const std::list<CLONES::Mutants::CellId> &cell_ids) const
    {
        return ForestCore::get_coalescent_cells(
            static_cast<const CLONES::Mutations::PhylogeneticForest &>(*this), cell_ids);
    }

    PhylogeneticForest
    get_subforest_for(const std::vector<std::string> &sample_names) const;

    Rcpp::DataFrame get_absolute_chromosome_positions() const;

    Rcpp::DataFrame get_germline_SIDs() const;

    inline const std::map<CLONES::Mutations::SID, std::string>& get_driver_codes() const
    {
        return driver_codes;
    }

    inline Rcpp::DataFrame get_sampled_cell_SIDs() const
    {
        return get_sampled_cell_SIDs(false);
    }

    Rcpp::DataFrame get_sampled_cell_SIDs(const bool with_germline) const;

    Rcpp::DataFrame get_sampled_cell_CNAs() const;

    Rcpp::List get_first_occurrence(const SEXP &mutation) const;

    Rcpp::List get_timed_exposures() const;

    Rcpp::List get_bulk_allelic_fragmentation() const;

    Rcpp::List get_bulk_allelic_fragmentation(const std::string &sample_name) const;

    Rcpp::List get_cell_allelic_fragmentation() const;

    inline std::string get_reference_path() const
    {
        return to_string(reference_path);
    }

    void set_reference_path(const std::string& reference_path);

    inline const GermlineSubject &get_germline_subject() const
    {
        return germline_subject;
    }

    inline Rcpp::List get_germline_subject_df() const
    {
        return germline_subject.get_dataframe();
    }

    inline void partition_samples(const SEXP &labelling_function)
    {
        ForestCore::partition_samples_in_forest(*this, labelling_function);
    }

    inline void save(const std::string &filename) const
    {
        save(filename, false);
    }

    void save(const std::string &filename, const bool quiet) const;

    inline static PhylogeneticForest load(const std::string &filename)
    {
        return load(filename, false);
    }

    static PhylogeneticForest load(const std::string &filename, const bool quiet);

    void show() const;

    friend class MutationEngine;
};

RCPP_EXPOSED_CLASS(PhylogeneticForest)
RCPP_EXPOSED_CLASS_NODECL(PhylogeneticForest::const_node)
using PhylogeneticForestNode = PhylogeneticForest::const_node;


#endif // __PROCESS_PHYLOGENETIC_FOREST__
