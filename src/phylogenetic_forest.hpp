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
#include "sampled_cell.hpp"

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

        for (const auto [mutation, cell_ids]: first_occurrences) {
            for (const auto cell_id: cell_ids) {
                ++(new_mutations[cell_id]);
            }
        }

        return new_mutations;
    }

    std::map<CLONES::Mutants::CellId, size_t>
    get_total_mutations(const std::map<CLONES::Mutants::CellId, size_t>& new_mutations) const;

  public:
    using base_type = CLONES::Mutations::PhylogeneticForest;

    class const_node
    {
    public:
        const_node(const PhylogeneticForest* forest, const base_type::const_node& node):
            node{node}, forest{forest}
        {}

        const_node(const PhylogeneticForest& forest, const base_type::const_node& node):
            const_node{&forest, node}
        {}

        const_node(const PhylogeneticForest* forest, const CLONES::Mutants::CellId cell_id):
            node{forest, cell_id}, forest{forest}
        {}

        const_node(const PhylogeneticForest& forest, const CLONES::Mutants::CellId cell_id):
            const_node{&forest, cell_id}
        {}

        inline const_node parent() const
        {
            return {forest, node.parent()};
        }

        std::vector<const_node> children() const
        {
            std::vector<const_node> children;

            for (const auto child: node.children()) {
                children.emplace_back(forest, child);
            }

            return children;
        }

        inline CLONES::Mutants::CellId get_id() const
        {
            return node.get_id();
        }

        inline size_t get_height() const
        {
            return node.height();
        }

        inline bool is_leaf() const
        {
            return node.is_leaf();
        }

        inline bool is_root() const
        {
            return node.is_root();
        }

        inline CLONES::Time get_birth_time() const
        {
            return node.get_birth_time();
        }

        inline CLONES::Time get_death_time() const
        {
            return node.get_death_time();
        }

        inline CLONES::Time get_life_span() const
        {
            return get_death_time()-get_birth_time();
        }

        inline CLONES::Mutants::SpeciesId get_species_id() const
        {
            return node.get_species_id();
        }

        inline CLONES::Mutants::MutantId get_mutant_id() const
        {
            return node.get_mutant_id();
        }

        inline std::string get_species_name() const
        {
            return node.get_species_name();
        }

        inline std::string get_mutant_name() const
        {
            return node.get_mutant_name();
        }

        Rcpp::List arising_mutations() const;
    private:
        base_type::const_node node;
        const PhylogeneticForest* forest;
    };

    PhylogeneticForest();

    std::vector<const_node> get_roots() const;

    inline Rcpp::List get_nodes() const
    {
        return ForestCore::get_nodes(
            static_cast<const CLONES::Mutations::PhylogeneticForest &>(*this));
    }

    Rcpp::List get_mutation_statistics() const;

    const std::list<CLONES::Mutants::CellId> &
    get_cell_ids_in(const std::string &sample_name) const;

    Rcpp::List get_samples_info() const;

    Rcpp::List get_driver_mutations() const;

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

    Rcpp::List get_absolute_chromosome_positions() const;

    Rcpp::List get_germline_SIDs() const;

    Rcpp::List get_sampled_cell_SIDs() const;

    Rcpp::List get_cell_SIDs(const CLONES::Mutants::CellId &cell_id) const;

    Rcpp::List get_sampled_cell_CNAs() const;

    Rcpp::List get_cell_CNAs(const CLONES::Mutants::CellId &cell_id) const;

    Rcpp::List get_first_occurrence(const SEXP &mutation) const;

    Rcpp::List get_timed_exposures() const;

    Rcpp::List get_bulk_allelic_fragmentation() const;

    Rcpp::List get_bulk_allelic_fragmentation(const std::string &sample_name) const;

    Rcpp::List get_cell_allelic_fragmentation() const;

    inline std::filesystem::path get_reference_path() const { return reference_path; }

    void set_reference_path(const std::filesystem::path reference_path);

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

    inline void save(const std::string &filename) const { save(filename, false); }

    void save(const std::string &filename, const bool quiet) const;

    inline static PhylogeneticForest load(const std::string &filename)
    {
        return load(filename, false);
    }

    static PhylogeneticForest load(const std::string &filename, const bool quiet);

    static Rcpp::List
    get_SID_dataframe(const CLONES::Mutations::CellGenomeMutations &cell_mutations);

    void show() const;

    friend class MutationEngine;
    friend class PhylogeneticForest::const_node;
};

RCPP_EXPOSED_CLASS(PhylogeneticForest)
RCPP_EXPOSED_CLASS_NODECL(PhylogeneticForest::const_node)
using PhylogeneticForestNode = PhylogeneticForest::const_node;


#endif // __PROCESS_PHYLOGENETIC_FOREST__
