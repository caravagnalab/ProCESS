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
    using CLONES_Tour = typename CLONES::Mutations::LabelTour<FOREST, LabellingFunctor<FOREST>>;

    std::shared_ptr<CLONES_Tour> tour;
    std::shared_ptr<typename CLONES_Tour::const_iterator> tour_it;

public:
    using functor_type = LabellingFunctor<FOREST>;
    using label_type = functor_type::label_type;

    LabelTour(const FOREST& forest, const functor_type& functor,
              const label_type& init_label, const bool only_leaves):
        tour{std::make_shared<CLONES_Tour>(forest, functor, init_label, only_leaves)},
        tour_it{std::make_shared<typename CLONES_Tour::const_iterator>(tour->begin())}
    {}

    Rcpp::List get_value() const
    {
        using namespace Rcpp;
        const auto& value = *(*tour_it);

        return DataFrame::create(_["cell_id"] = value.first,
                                 _["label"] = value.second);
    }

    inline void step()
    {
        tour_it->operator++();
    }

    inline bool done() const
    {
        return tour_it->is_end();
    }
};

template<typename FOREST>
LabelTour<FOREST> get_label_tour(const FOREST& forest,
                                 const Rcpp::Function& R_functor,
                                 const Rcpp::RObject& init_label,
                                 const bool only_leaves)
{
    LabellingFunctor<FOREST> C_functor{R_functor};

    return {forest, C_functor, init_label, only_leaves};
}

/*
Rcpp::RObject get_label_tour(const Rcpp::RObject& forest,
                             const Rcpp::Function& R_functor,
                             const Rcpp::RObject& init_label,
                             const bool only_leaves);
*/

#endif // __PROCESS_FOREST_LABELLING__
