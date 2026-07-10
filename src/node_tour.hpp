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

#ifndef __PROCESS_NODE_TOUR__
#define __PROCESS_NODE_TOUR__

#include <Rcpp.h>

#include <label_tour.hpp>

#include "cell_mutations.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"

#define REGISTER_TOUR_COMMON_FIELDS(ClassType)                  \
        .property("node", &ClassType::get_node,                 \
            "The node pointed by the iterator")                 \
        .property("label", &ClassType::get_label,               \
            "The label of the node pointed by the iterator")    \
        .method("step", &ClassType::step,                       \
            "Go to the next node in the tour")                  \
        .property("done", &ClassType::done,                     \
            "Test whether the tour ended")                      \
        .method("show", &ClassType::show,                       \
            "Show the ClassType object")                        \


template<typename FOREST>
class LabellingFunctor
{
    const Rcpp::Function R_function;

public:
    using label_type = Rcpp::RObject;
    using const_node = typename FOREST::const_node;

    LabellingFunctor(const Rcpp::Function& R_function):
        R_function{R_function}
    {}

    inline label_type operator()(const label_type& parent_label,
                                 const const_node& node) const
    {
        return R_function(parent_label, node);
    }
};

template<typename FOREST>
class NodeFunctor
{
public:
    using const_node = typename FOREST::const_node;
    using label_type = const_node;

    NodeFunctor()
    {}

    inline label_type operator()(const label_type& parent_label,
                                 const const_node& node) const
    {
        return node;
    }
};

template<typename FOREST>
    requires CLONES::Mutations::isForest<FOREST>
class NodeTour
{
public:
    using node_functor_type = NodeFunctor<FOREST>;
    using genome_functor_type =  CLONES::Mutations::MutationLabellingFunctor<GenomeMutations>;
    using labelling_functor_type = LabellingFunctor<FOREST>;

    using node_tour_type = typename CLONES::Mutations::LabelTour<FOREST, node_functor_type>;
    using genome_tour_type = typename CLONES::Mutations::LabelTour<CLONES::Mutations::PhylogeneticForest, genome_functor_type>;
    using label_tour_type = typename CLONES::Mutations::LabelTour<FOREST, labelling_functor_type>;

    using label_type = typename labelling_functor_type::label_type;

    NodeTour(const FOREST& forest, const labelling_functor_type& label_functor,
              const label_type& init_label, const bool only_leaves = false,
              const bool with_genome = false):
        NodeTour(forest, only_leaves, with_genome)
    {
        label_tour = std::make_shared<label_tour_type>(forest, label_functor, init_label, only_leaves);

        reset_tour(label_tour, label_tour_it);
    }

    NodeTour(const FOREST& forest, const bool only_leaves = false,
              const bool with_genome = false):
        label_tour{nullptr}
    {
        if (with_genome) {
            if constexpr (std::is_base_of_v<PhylogeneticForest, FOREST>) {
                genome_functor_type g_functor(true);

                if (!std::filesystem::exists(forest.get_reference_path())) {
                    Rcpp::stop("The forest reference FASTA file \""
                                + to_string(forest.get_reference_path())
                                + "\" does not exists. Reset it in forest by "
                                + "using `PhylogeneticForest::set_reference_path()`.");
                }

                GenomeMutations init_mutations{forest.get_reference_path(),
                                               forest.get_germline_mutations()};

                genome_tour = std::make_shared<genome_tour_type>(forest, g_functor,
                                                                 init_mutations,
                                                                 only_leaves);

                reset_tour(genome_tour, genome_tour_it);
            } else {
                Rcpp::stop("Only `PhylogeneticForest` objects supports genome generation.");
            }
        }

        typename FOREST::const_node node;

        NodeFunctor<FOREST> n_functor;

        node_tour = std::make_shared<node_tour_type>(forest, n_functor, node, only_leaves);

        reset_tour(node_tour, node_tour_it);
    }

    NodeTour(const NodeTour<FOREST>& orig):
        node_tour{orig.node_tour}, genome_tour{orig.genome_tour}, label_tour{orig.label_tour},
        node_tour_it{copy_iterator(orig.node_tour_it)},
        genome_tour_it{copy_iterator(orig.genome_tour_it)},
        label_tour_it{copy_iterator(orig.label_tour_it)}
    {}

    inline const FOREST::const_node& get_node() const
    {
        return get_value(node_tour_it, "node", "");
    }

    inline const label_type& get_label() const
    {
        return get_value(label_tour_it, "label", "labelling_functor");
    }

    inline const GenomeMutations& get_genome() const
    {
        return get_value(genome_tour_it, "genome", "with_genomes");
    }

    void step()
    {
        if (node_tour_it != nullptr) {
            node_tour_it->operator++();
        }

        if (label_tour_it != nullptr) {
            label_tour_it->operator++();
        }
    
        if (genome_tour_it != nullptr) {
            genome_tour_it->operator++();
        }
    }

    inline bool done() const
    {
        if (node_tour_it != nullptr) {
            return node_tour_it->is_end();
        }

        return true;
    }

    void reset_tour()
    {
        reset_tour(node_tour, node_tour_it);
        reset_tour(genome_tour, genome_tour_it);
        reset_tour(label_tour, label_tour_it);
    }

    void show() const
    {
        Rcpp::Rcout << get_demangled_type_name(typeid(FOREST)) 
                    << "NodeTour(with_labels = "
                    << (label_tour != nullptr?"TRUE":"FALSE");

        if constexpr (std::is_base_of_v<PhylogeneticForest, FOREST>) {
            Rcpp::Rcout << ", with_genomes = "
                        << (genome_tour != nullptr?"TRUE":"FALSE");
        }

        Rcpp::Rcout << ")" << std::endl;
    }
private:
    std::shared_ptr<node_tour_type> node_tour;
    std::shared_ptr<genome_tour_type> genome_tour;
    std::shared_ptr<label_tour_type> label_tour;

    std::unique_ptr<typename node_tour_type::const_iterator> node_tour_it;
    std::unique_ptr<typename genome_tour_type::const_iterator> genome_tour_it;
    std::unique_ptr<typename label_tour_type::const_iterator> label_tour_it;

    template<typename LABEL_TOUR>
    static void reset_tour(const std::shared_ptr<LABEL_TOUR>& tour,
                           std::unique_ptr<typename LABEL_TOUR::const_iterator>& tour_it)
    {
        if (tour != nullptr) {
            tour_it = std::make_unique<typename LABEL_TOUR::const_iterator>(tour->begin());
        } else {
            tour_it = nullptr;
        }
    }
protected:
    NodeTour():
        node_tour{nullptr}, genome_tour{nullptr}, label_tour{nullptr}
    {
        reset_tour();
    }

    template<typename LABEL_ITERATOR>
    static const typename LABEL_ITERATOR::label_type&
    get_value(const std::unique_ptr<LABEL_ITERATOR>& iterator,
              const char* name, const char* parameter_name)
    {
        if (iterator == nullptr) {
            std::string str_name{name};
            if constexpr (!std::is_base_of_v<PhylogeneticForest, FOREST>) {
                Rcpp::stop("The field `" + str_name
                           + "` exclusively available in `PhylogeneticForest` "
                           + "node tours.");
            } else {
                Rcpp::stop("The field `" + str_name
                           + "` not available. Use the `get_node_tour()`'s "
                           + "parameter \"" + parameter_name
                           + "\" to enable it.");
            }
        }

        return (*(*iterator)).second;
    }

    template<typename LABEL_ITERATOR>
    static std::unique_ptr<LABEL_ITERATOR> copy_iterator(const std::unique_ptr<LABEL_ITERATOR>& iterator)
    {
        if (iterator == nullptr) {
            return nullptr;
        }

        return std::make_unique<LABEL_ITERATOR>(*iterator);
    }
};

SEXP get_node_tour(const SEXP& forest, const Rcpp::RObject& R_functor,
                   const Rcpp::RObject& init_label, const bool only_leaves,
                   const bool with_genomes);


RCPP_EXPOSED_CLASS_NODECL(NodeTour<SampleForest>)
RCPP_EXPOSED_CLASS_NODECL(NodeTour<PhylogeneticForest>)

using SampleForestNodeTour = NodeTour<SampleForest>;
using PhylogeneticForestNodeTour = NodeTour<PhylogeneticForest>;

#endif // __PROCESS_NODE_TOUR__
