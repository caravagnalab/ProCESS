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

#ifndef __PROCESS_SAMPLES_FOREST__
#define __PROCESS_SAMPLES_FOREST__

#include <vector>

#include <Rcpp.h>

#include <phylogenetic_forest.hpp>

#include "forest.hpp"
#include "forest_labelling.hpp"

class MutationEngine;

class SampleForest : public CLONES::Mutants::DescendantForest
{
    SampleForest();

    Rcpp::List get_nodes(const std::vector<CLONES::Mutants::CellId> &cell_ids) const;

  public:
    using base_type = CLONES::Mutants::DescendantForest;

    class const_node
    {
    public:
        const_node(const base_type::const_node& node):
            node{node}
        {}

        const_node(const SampleForest* forest, const CLONES::Mutants::CellId cell_id):
            node{forest, cell_id}
        {}

        const_node(const SampleForest& forest, const CLONES::Mutants::CellId cell_id):
            const_node{&forest, cell_id}
        {}

        inline const_node parent() const
        {
            return {node.parent()};
        }

        std::vector<const_node> children() const
        {
            std::vector<const_node> children;

            for (const auto child: node.children()) {
                children.emplace_back(child);
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
    private:
        base_type::const_node node;
    };

    SampleForest(const CLONES::Mutants::Evolutions::TissueSimulation &simulation);

    std::vector<const_node> get_roots() const;

    inline Rcpp::List get_nodes() const
    {
        return ForestCore::get_nodes(
            static_cast<const CLONES::Mutants::DescendantForest &>(*this));
    }

    Rcpp::List get_samples_info() const;

    inline Rcpp::List get_species_info() const
    {
        return ForestCore::get_species_info(
            static_cast<const CLONES::Mutants::DescendantForest &>(*this));
    }

    inline Rcpp::List get_coalescent_cells() const
    {
        return ForestCore::get_coalescent_cells(
            static_cast<const CLONES::Mutants::DescendantForest &>(*this));
    }

    inline Rcpp::List
    get_coalescent_cells(const std::list<CLONES::Mutants::CellId> &cell_ids) const
    {
        return ForestCore::get_coalescent_cells(
            static_cast<const CLONES::Mutants::DescendantForest &>(*this), cell_ids);
    }

    inline std::list<std::list<CLONES::Mutants::CellId>> get_sticks() const
    {
        return CLONES::Mutants::DescendantForest::get_sticks();
    }

    inline std::list<std::list<CLONES::Mutants::CellId>>
    get_sticks(const double birth_threshold) const
    {
        return CLONES::Mutants::DescendantForest::get_sticks(birth_threshold);
    }

    SampleForest get_subforest_for(const std::vector<std::string> &sample_names) const;

    void save(const std::string &filename) const;

    static SampleForest load(const std::string &filename);

    void show() const;

    friend class MutationEngine;
};

inline SampleForest load_samples_forest(const std::string &filename)
{
    Rcpp::warning("`load_samples_forest()` is deprecated. "
                  "Please use `load_sample_forest()` instead.");

    return SampleForest::load(filename);
}

RCPP_EXPOSED_CLASS(SampleForest)
RCPP_EXPOSED_CLASS_NODECL(SampleForest::const_node)
RCPP_EXPOSED_CLASS_NODECL(LabelTour<SampleForest>)

using SampleForestNode = SampleForest::const_node;
using SampleForestLabelTour = LabelTour<SampleForest>;

#endif // __PROCESS_SAMPLES_FOREST__
