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

#ifndef __PROCESS_FOREST__
#define __PROCESS_FOREST__

#include <sstream>
#include <vector>

#include <Rcpp.h>

#include <cell.hpp>
#include <label_tour.hpp>
#include <mutant_properties.hpp>
#include <mutation_list.hpp>

#include "cell_mutations.hpp"
#include "utility.hpp"

#define REGISTER_FOREST_COMMON_FIELD(ClassType)                                     \
        .method("represents_cell", &ClassType::represents_cell,                     \
                "Test whether the forest has a node representing a cell")           \
        .method("get_nodes",                                                        \
                (Rcpp::List (ClassType::*)() const)(&ClassType::get_nodes),         \
                "Get the nodes of the forest")                                      \
        .method("get_node",                                                         \
                 (ClassType::const_node (ClassType::*)(const SEXP) const)(          \
                    &ClassType::get_node),                                          \
                "Get node corresponding to a cell")                                 \
        .method("get_coalescent_cells",                                             \
                (List (ClassType::*)(const std::list<CLONES::Mutants::CellId> &)    \
                     const)(&ClassType::get_coalescent_cells),                      \
                "Get the most recent common ancestor of some cells")                \
        .method("get_coalescent_cells",                                             \
                (List (ClassType::*)() const)(&ClassType::get_coalescent_cells),    \
                "Get the most recent common ancestor of all forest leaves")         \
        .method("get_subforest_for", &ClassType::get_subforest_for,                 \
                "Get the sub-forest for some of the original samples")              \
        .method("get_samples_info", &ClassType::get_samples_info,                   \
                "Get some pieces of information about the samples")                 \
        .method("get_species_info", &ClassType::get_species_info,                   \
                "Get the recorded species")                                         \
        .method("get_sticks",                                                       \
                (std::list<std::list<CLONES::Mutants::CellId>> (ClassType::*)(      \
                    const double) const)(&ClassType::get_sticks),                   \
                "Get the forest sticks")                                            \
        .method("get_sticks",                                                       \
                (std::list<std::list<CLONES::Mutants::CellId>> (ClassType::*)()     \
                     const)(&ClassType::get_sticks),                                \
                "Get the forest sticks")                                            \
        .method("save",                                                             \
                (void (ClassType::*)(const std::string &, const bool) const)        \
                    &ClassType::save, "Save a sample forest")                       \
        .method("save",                                                             \
                (void (ClassType::*)(const std::string &) const)                    \
                    &ClassType::save, "Save a sample forest")                       \
        .method("show", &ClassType::show, "Describe the forest")                    \


#define REGISTER_NODE_COMMON_FIELDS(ClassType)                      \
        .property("cell_id", &ClassType::get_id,                    \
            "The identifier of the cell associated to the "         \
            "node (Read-only)")                                     \
        .property("parent", &ClassType::parent,                     \
            "The node's parent (Read-only)")                        \
        .property("children", &ClassType::children,                 \
            "The node's children (Read-only)")                      \
        .property("is_root", &ClassType::is_root,                   \
            "A Boolean flag that is set to `TRUE` iff the "         \
            "node is a root (Read-only)")                           \
        .property("is_leaf", &ClassType::is_leaf,                   \
            "A Boolean flag that is set to TRUE iff the "           \
            "node is a leaf (Read-only)")                           \
        .property("sample_name", &ClassType::get_sample_name,       \
            "The name of the sample, when available, that "         \
            "collected the cell associated to the node "            \
            "(Read-only)")                                          \
        .property("birth_time", &ClassType::get_birth_time,         \
            "The birth time of the cell associated to the "         \
            "node (Read-only)")                                     \
        .property("death_time", &ClassType::get_death_time,         \
            "The death time of the cell associated to the "         \
            "node (Read-only)")                                     \
        .property("life_span", &ClassType::get_life_span,           \
            "The life span of the cell associated to the "          \
            "node (Read-only)")                                     \
        .property("species_id", &ClassType::get_species_id,         \
            "The identifier of the species of the cell "            \
            "associated to the node (Read-only)")                   \
        .property("species_name", &ClassType::get_species_name,     \
            "The name of the species of the cell associated"        \
            "to the node (Read-only)")                              \
        .property("epistate_name", &ClassType::get_epistate_name,   \
            "The name of the epistate of the cell associated "      \
            "to the node(Read-only)")                               \
        .property("mutant_id", &ClassType::get_mutant_id,           \
            "The identifier of the mutant of the cell "             \
            "associated to the node (Read-only)")                   \
        .property("mutant_name", &ClassType::get_mutant_name,       \
            "The name of the mutant of the cell associated "        \
            "to the node (Read-only)")                              \
        .method("show", &ClassType::show, "Show the ClassType")     \


template <typename T>
concept hasArisingMutationMethod = requires(const T instance) {
    { instance.arising_mutations() } -> std::same_as<CLONES::Mutations::MutationList>;
};

class ForestCore
{
    template<typename FOREST>
    class LabellingFunction
    {
        FOREST const& forest;
        Rcpp::Function const& R_function;

    public:
        LabellingFunction(const FOREST& forest,
                          const Rcpp::Function& R_function):
            forest{forest}, R_function{R_function}
        {}

        std::string operator()(const CLONES::Mutants::CellId& cell_id) const
        {
            auto node = Rcpp::wrap(typename FOREST::const_node{forest, cell_id});

            return Rcpp::as<std::string>(R_function(node));
        }
    };
public:
    inline static void fill_mutation_list(Rcpp::DataFrame &df, const CLONES::Mutations::MutationList& mutations,
                                          const std::string& mutant_name,
                                          const std::map<CLONES::Mutations::SID, std::string>& driver_codes)
    {
        size_t i{0};

        fill_mutation_list(df, mutations, mutant_name, driver_codes, i);
    }

    inline static void fill_mutation_list(Rcpp::DataFrame &df, const CLONES::Mutations::MutationList& mutations,
                                          const std::string& mutant_name,
                                          const std::map<CLONES::Mutations::SID, std::string>& driver_codes,
                                          size_t &i)
    {
        size_t mut_i{1};

        fill_mutation_list(df, mutations, mutant_name, driver_codes, i, mut_i);
    }

    static void fill_mutation_list(Rcpp::DataFrame &df, const CLONES::Mutations::MutationList& mutations,
                                   const std::string& mutant_name,
                                   const std::map<CLONES::Mutations::SID, std::string>& driver_codes,
                                   size_t &i, size_t &mut_i);

    template<typename FOREST>
#if defined(__clang__)
      /* clang and GCC have different ways of handling concept
         verification on incomplete types. Clang uses lazy evaluation
         and let us requiring a concept for a template parameter even
         if the parameter is not completely defined. GCC, instead,
         requires a fully defined parameter to verify the concept.
         To defer evaluation of concepts, we can move concept
         requirements from template parameter to constructors.
         However, this solution does not work in clang which,
         apparently, defers evaluation on template parameters, but
         not over the constructors. At the moment, there are two
         solutions: we can either 1) distinguish the two compilers
         and apply ad-hoc code or 2) remove the isForest
         requirement tout court */

      requires CLONES::Mutations::isForest<FOREST>
#endif
    class const_node
    {
        constexpr inline static std::string get_forest_class_name()
        {
            return get_demangled_type_name(typeid(FOREST));
        }
    public:
        const_node()
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : node{}, forest{nullptr}
        {}

        const_node(const FOREST* forest, const FOREST::base_type::const_node& node)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : node{node}, forest{forest}
        {}

        const_node(const FOREST& forest, const FOREST::base_type::const_node& node)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : const_node<FOREST>{&forest, node}
        {}

        const_node(const FOREST* forest, const CLONES::Mutants::CellId cell_id)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : node{forest, cell_id}, forest{forest}
        {}

        const_node(const FOREST& forest, const CLONES::Mutants::CellId cell_id)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : const_node<FOREST>{&forest, cell_id}
        {}

        const_node(const FOREST* forest, const SEXP cell_id)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : forest{forest}
        {
            using CellId = CLONES::Mutants::CellId;

            const auto C_cell_id = FromSEXP::get<CellId>(cell_id, "parameter `cell_id`",
                                                         "an integer value");

            node = typename FOREST::base_type::const_node(forest, C_cell_id);
        }

        const_node(const FOREST& forest, const SEXP cell_id)
#if !defined(__clang__) && defined(__GNUC__)
            requires CLONES::Mutations::isForest<FOREST>
#endif
            : const_node<FOREST>{&forest, cell_id}
        {}

        const_node<FOREST> parent() const
        {
            if (is_root()) {
                Rcpp::stop("The node is a forest root and has no parent.");
            }
            return {forest, node.parent()};
        }

        std::vector<const_node<FOREST>> children() const
        {
            std::vector<const_node<FOREST>> children;

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

        inline std::string get_epistate_name() const
        {
            return node.get_epistate_name();
        }

        Rcpp::String get_sample_name() const
        {
            using namespace Rcpp;

            if (is_leaf()) {
                return node.get_sample().get_name();
            }

            return NA_STRING;
        }

        void show() const
        {
            Rcpp::Rcout << get_forest_class_name() << "Node(";
            if (forest != nullptr) {
                Rcpp::Rcout << "cell_id = "
                            << static_cast<size_t>(get_id())
                            << ", species = \"" << get_species_name()
                            << "\"";
            }

            Rcpp::Rcout << ")" << std::endl;
        }

        GenomeMutations get_genome() const
        {
            node.assert_initialized();

            const auto somatic = node.cell_mutations(true, false);

            return {forest->get_reference_path(),
                    forest->get_germline_mutations(),
                    somatic};
        }

        Rcpp::List arising_mutations() const
        {
            using namespace Rcpp;

            const std::string mutant_name{get_mutant_name()};

            if constexpr (hasArisingMutationMethod<typename FOREST::base_type::const_node>) {
                const auto& mutations = node.arising_mutations();

                const size_t num_of_rows{mutations.size()};

                CharacterVector chrs(num_of_rows), types(num_of_rows),
                    CNA_types(num_of_rows), refs(num_of_rows), alts(num_of_rows), codes(num_of_rows),
                    natures(num_of_rows), causes(num_of_rows);
                IntegerVector starts(num_of_rows), ends(num_of_rows), application_order(num_of_rows),
                    alleles(num_of_rows), src_alleles(num_of_rows);

                DataFrame df = DataFrame::create(_["order"] = application_order,
                                                 _["type"] = types, _["CNA_type"] = CNA_types,
                                                 _["chr"] = chrs, _["start"] = starts,
                                                 _["end"] = ends, _["ref"] = refs,
                                                 _["alt"] = alts, _["allele"] = alleles,
                                                 _["src.allele"] = src_alleles,
                                                 _["nature"] = natures, _["cause"] = causes,
                                                 _["code"] = codes);

                fill_mutation_list(df, mutations, mutant_name, forest->get_driver_codes());

                return df;
            }

            std::ostringstream oss;
            oss << get_demangled_type_name(typeid(FOREST)) << "Node does not "
                << "support the method `arising_mutations()`.";

            stop(oss.str());
        }
    private:
        FOREST::base_type::const_node node;
        FOREST const* forest;
    };

    inline static bool represents_cell(const CLONES::Mutants::DescendantForest& forest, const SEXP cell_id)
    {
        using CellId = CLONES::Mutants::CellId;

        const auto C_cell_id = FromSEXP::get<CellId>(cell_id, "parameter `cell_id`",
                                                     "an integer value");

        const auto& cells = forest.get_cells();

        return (cells.find(C_cell_id) != cells.end());
    }

    template <typename CPP_FOREST>
    static Rcpp::List get_nodes(const CPP_FOREST &forest,
                                const std::vector<CLONES::Mutants::CellId> &cell_ids)
    {
        using namespace Rcpp;
        using namespace CLONES::Mutants;

        IntegerVector ids(cell_ids.size()), ancestors(cell_ids.size()),
                      depths(cell_ids.size());
        CharacterVector mutants(cell_ids.size()), epi_states(cell_ids.size()),
            sample_names(cell_ids.size());
        NumericVector birth(cell_ids.size());

        const auto node_depths = forest.get_node_depths();

        bool with_epigenetics{false};
        size_t i{0};
        for (const auto &cell_id : cell_ids) {
            ids[i] = cell_id;
            auto cell_node = forest.get_node(cell_id);
            if (cell_node.is_root()) {
                ancestors[i] = NA_INTEGER;
            } else {
                ancestors[i] = cell_node.parent().get_id();
            }

            depths[i] = node_depths.at(cell_id);
            mutants[i] = cell_node.get_mutant_name();
            const auto epistate = cell_node.get_epistate_name();
            if (epistate != "") {
                epi_states[i] = epistate;
                with_epigenetics = true;
            }

            if (cell_node.is_leaf()) {
                sample_names[i] = cell_node.get_sample().get_name();
            } else {
                sample_names[i] = NA_STRING;
            }
            birth[i] = static_cast<const Cell &>(cell_node).get_birth_time();

            ++i;
        }

        if (with_epigenetics) {
            return DataFrame::create(_["cell_id"] = ids, _["ancestor"] = ancestors,
                                    _["depth"] = depths, _["mutant"] = mutants,
                                    _["epistate"] = epi_states, _["sample"] = sample_names,
                                    _["birth_time"] = birth);
        }

        return DataFrame::create(_["cell_id"] = ids, _["ancestor"] = ancestors,
                                 _["depth"] = depths, _["mutant"] = mutants,
                                 _["sample"] = sample_names, _["birth_time"] = birth);
    }

    template <typename CPP_FOREST> static Rcpp::List get_nodes(const CPP_FOREST &forest)
    {
        std::vector<CLONES::Mutants::CellId> cell_ids;
        cell_ids.reserve(forest.num_of_nodes());
        for (const auto &[cell_id, cell] : forest.get_cells()) {
            cell_ids.push_back(cell_id);
        }

        return ForestCore::get_nodes<CPP_FOREST>(forest, cell_ids);
    }

    template <typename CPP_FOREST>
    static Rcpp::List get_species_info(const CPP_FOREST &forest)
    {
        using namespace Rcpp;

        size_t num_of_rows = forest.get_species_data().size();

        CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);

        using namespace CLONES::Mutants;

        bool with_epigenetics{false};
        size_t i{0};
        for (const auto &[species_id, species_data] : forest.get_species_data()) {
            mutant_names[i] = species_data.get_mutant_name();

            const auto epistate = species_data.get_epistate_name();
            if (epistate != "") {
                epi_states[i] = epistate;
                with_epigenetics = true;
            }

            ++i;
        }

        if (with_epigenetics) {
            return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states);
        }

        return DataFrame::create(_["mutant"] = mutant_names);
    }

    template <typename CPP_FOREST>
    static Rcpp::List get_coalescent_cells(const CPP_FOREST &forest)
    {
        auto coalencent_ids = forest.get_coalescent_cells();

        return ForestCore::get_nodes<CPP_FOREST>(forest, coalencent_ids);
    }

    template <typename CPP_FOREST>
    static Rcpp::List
    get_coalescent_cells(const CPP_FOREST &forest,
                         const std::list<CLONES::Mutants::CellId> &cell_ids)
    {
        auto coalencent_ids = forest.get_coalescent_cells(cell_ids);

        return ForestCore::get_nodes<CPP_FOREST>(forest, coalencent_ids);
    }

    template<typename FOREST>
    static void
    partition_samples_in_forest(FOREST& forest, const SEXP &labelling_function)
    {
        switch (TYPEOF(labelling_function)) {
        case CLOSXP:
        {
            Rcpp::Function l_function = Rcpp::as<Rcpp::Function>(labelling_function);

            LabellingFunction l_func{forest, l_function};

            static_cast<FOREST::base_type&>(forest).partition_samples(l_func);

            break;
        }
        default:
            Rcpp::stop("The parameter `labelling_function` must be a function.");
        }
    }

    static SEXP load_forest(const std::string &filename, const bool quiet);
};

#endif // __PROCESS_FOREST__
