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

class MutationEngine;

class SampleForest : public CLONES::Mutants::DescendantForest
{
    SampleForest();

    Rcpp::List get_nodes(const std::vector<CLONES::Mutants::CellId> &cell_ids) const;

  public:
    inline static const std::string file_format_header = "ProCESS Sample Forest";
    inline static const uint8_t file_format_number = 0;

    using base_type = CLONES::Mutants::DescendantForest;

    using const_node = ForestCore::const_node<SampleForest>;

    SampleForest(const CLONES::Mutants::Evolutions::TissueSimulation &simulation);

    std::vector<const_node> get_roots() const;

    inline Rcpp::List get_nodes() const
    {
        return ForestCore::get_nodes(
            static_cast<const CLONES::Mutants::DescendantForest &>(*this));
    }

    inline const_node get_node(const CLONES::Mutants::CellId &cell_id) const
    {
        return SampleForest::const_node(*this, cell_id);
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

    inline void save(const std::string &filename) const
    {
        save(filename, false);
    }

    void save(const std::string &filename, const bool quiet) const;

    static SampleForest load(const std::string &filename, const bool quiet);

    void show() const;

    friend class MutationEngine;
};

RCPP_EXPOSED_CLASS(SampleForest)
RCPP_EXPOSED_CLASS_NODECL(SampleForest::const_node)

using SampleForestNode = SampleForest::const_node;


#endif // __PROCESS_SAMPLES_FOREST__
