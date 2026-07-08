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

#ifndef __PROCESS_FOREST_LABELLING__
#define __PROCESS_FOREST_LABELLING__

#include <Rcpp.h>

#include <label_tour.hpp>

#include "cell_mutations.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"

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
    requires CLONES::Mutations::isForest<FOREST>
class LabelTour
{
public:
    using labelling_functor_type = LabellingFunctor<FOREST>;
    using genome_functor_type =  CLONES::Mutations::MutationLabellingFunctor<GenomeMutations>;

    using label_tour_type = typename CLONES::Mutations::LabelTour<FOREST, labelling_functor_type>;
    using genome_tour_type = typename CLONES::Mutations::LabelTour<CLONES::Mutations::PhylogeneticForest, genome_functor_type>;

    using label_type = typename labelling_functor_type::label_type;

    LabelTour(const FOREST& forest, const labelling_functor_type& label_functor,
              const label_type& init_label, const bool only_leaves = false,
              const bool with_genome = false):
        label_tour{std::make_shared<label_tour_type>(forest, label_functor, init_label, only_leaves)},
        genome_tour{nullptr}
    {
        if constexpr (std::is_base_of_v<PhylogeneticForest, FOREST>) {
            if (with_genome) {
                init_genome_tour(forest, only_leaves);
            }
        } else {
            if (with_genome) {
                Rcpp::stop("Only `PhylogeneticForest` objects supports genome generation.");
            }
        }

        reset_tour();
    }

    LabelTour(const PhylogeneticForest& forest, const bool only_leaves = false)
    {
        init_genome_tour(forest, only_leaves);

        reset_tour();
    }

    LabelTour(const LabelTour<FOREST>& orig):
        label_tour{orig.label_tour}, genome_tour{orig.genome_tour},
        label_tour_it{copy_iterator(orig.label_tour_it)},
        genome_tour_it{copy_iterator(orig.genome_tour_it)}
    {}

    Rcpp::List get_value() const
    {
        Rcpp::List list = Rcpp::List::create();

        set_in(list, label_tour_it, "label");
        set_in(list, genome_tour_it, "genome");

        return list;
    }

    void step()
    {
        if (label_tour_it != nullptr) {
            label_tour_it->operator++();
        }
    
        if (genome_tour_it != nullptr) {
            genome_tour_it->operator++();
        }
    }

    bool done() const
    {
        if (label_tour_it != nullptr) {
            return label_tour_it->is_end();
        }
    
        if (genome_tour_it != nullptr) {
            return genome_tour_it->is_end();
        }

        return true;
    }

    void reset_tour()
    {
        reset_tour(label_tour, label_tour_it);
        reset_tour(genome_tour, genome_tour_it);
    }
private:
    std::shared_ptr<label_tour_type> label_tour;
    std::shared_ptr<genome_tour_type> genome_tour;

    std::unique_ptr<typename label_tour_type::const_iterator> label_tour_it;
    std::unique_ptr<typename genome_tour_type::const_iterator> genome_tour_it;

    void init_genome_tour(const PhylogeneticForest& forest, const bool only_leaves)
    {
        genome_functor_type g_functor(true);

        if (!std::filesystem::exists(forest.get_reference_path())) {
            Rcpp::stop("The forest reference FASTA file \""
                        + to_string(forest.get_reference_path())
                        + "\" does not exists. Reset it in forest by "
                        + "using `PhylogeneticForest::set_reference_path()`.");
        }

        GenomeMutations init_mutations{forest.get_reference_path(), forest.get_germline_mutations()};

        genome_tour = std::make_shared<genome_tour_type>(forest, g_functor, init_mutations, only_leaves);
    }

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
    LabelTour():
        label_tour{nullptr}, genome_tour{nullptr}
    {
        reset_tour();
    }

    void reset_tour(std::shared_ptr<label_tour_type> label_tour, std::shared_ptr<genome_tour_type> genome_tour)
    {
        this->label_tour = label_tour;
        this->genome_tour = genome_tour;

        reset_tour();
    }

    template<typename LABEL_ITERATOR>
    static void set_in(Rcpp::List& list, const std::unique_ptr<LABEL_ITERATOR>& iterator, const char* name_in_list)
    {  
        if (iterator != nullptr) {
            const auto& value = *(*iterator);
    
            list["cell_id"] = value.first;
            list[name_in_list] = value.second;
        }
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

SEXP get_label_tour(const SEXP& forest, const Rcpp::Function& R_functor,
                    const Rcpp::RObject& init_label, const bool only_leaves,
                    const bool with_genomes);

SEXP get_genome_tour(const SEXP& forest, const bool only_leaves);


RCPP_EXPOSED_CLASS_NODECL(LabelTour<SampleForest>)
RCPP_EXPOSED_CLASS_NODECL(LabelTour<PhylogeneticForest>)

using SampleForestLabelTour = LabelTour<SampleForest>;
using PhylogeneticForestLabelTour = LabelTour<PhylogeneticForest>;

#endif // __PROCESS_FOREST_LABELLING__
