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

#include "node_tour.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"

template<typename FOREST>
inline SEXP get_node_tour_template(const SEXP& forest,
                                   const Rcpp::RObject& R_functor,
                                   const Rcpp::RObject& init_label,
                                   const bool only_leaves,
                                   const bool with_genomes)
{
    using namespace Rcpp;

    if constexpr (!std::is_base_of_v<PhylogeneticForest, FOREST>) {
        if (with_genomes) {
            Rcpp::stop("Genome generation is exclusively supported by `PhylogeneticForest`.");
        }
    }

    S4 forest_s4obj(forest);

    Environment env(forest_s4obj);

    XPtr<FOREST> forest_ptr(env.get(".pointer"));

    if (R_functor == R_NilValue) {
        NodeTour<FOREST> label_tour{*forest_ptr, only_leaves,
                                    with_genomes};

        return wrap(label_tour);
    }

    LabellingFunctor<FOREST> C_functor{as<Function>(R_functor)};

    NodeTour<FOREST> label_tour{*forest_ptr, C_functor,
                                init_label, only_leaves,
                                with_genomes};

    return wrap(label_tour);
}


SEXP get_node_tour(const SEXP& forest, const Rcpp::RObject& R_functor,
                   const Rcpp::RObject& init_label, const bool only_leaves,
                   const bool with_genomes)
{
    using namespace Rcpp;

    switch (TYPEOF(forest)) {
        case S4SXP:
        {
            S4 s4obj(forest);

            if (s4obj.is("Rcpp_SampleForest")) {

                return get_node_tour_template<SampleForest>(s4obj, R_functor,
                                                            init_label, only_leaves,
                                                            with_genomes);
            }

            if (s4obj.is("Rcpp_PhylogeneticForest")) {

                return get_node_tour_template<PhylogeneticForest>(s4obj, R_functor,
                                                                  init_label, only_leaves,
                                                                  with_genomes);
            }
        }
        default:
            Rcpp::stop("Unsupported forest type: expected either a `SampleForest`"
                       " or a `PhylogeneticForest` object.");
    }
}
