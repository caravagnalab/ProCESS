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

#include "forest_labelling.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"

template<typename FOREST>
inline SEXP get_label_tour_template(const SEXP& forest,
                                    const Rcpp::Function& R_functor,
                                    const Rcpp::RObject& init_label,
                                    const bool only_leaves,
                                    const bool with_genomes)
{
    using namespace Rcpp;

    LabellingFunctor<FOREST> C_functor{R_functor};
    S4 forest_s4obj(forest);

    Environment env(forest_s4obj);

    XPtr<FOREST> forest_ptr(env.get(".pointer"));

    LabelTour<FOREST> label_tour{*forest_ptr, C_functor,
                                 init_label, only_leaves,
                                 with_genomes};

    return wrap(label_tour);
}


SEXP get_label_tour(const SEXP& forest, const Rcpp::Function& R_functor,
                    const Rcpp::RObject& init_label, const bool only_leaves,
                    const bool with_genomes)
{
    using namespace Rcpp;

    switch (TYPEOF(forest)) {
        case S4SXP:
        {
            S4 s4obj(forest);

            if (s4obj.is("Rcpp_SampleForest")) {

                return get_label_tour_template<SampleForest>(s4obj, R_functor,
                                                             init_label, only_leaves,
                                                             with_genomes);
            }

            if (s4obj.is("Rcpp_PhylogeneticForest")) {

                return get_label_tour_template<PhylogeneticForest>(s4obj, R_functor,
                                                                   init_label, only_leaves,
                                                                   with_genomes);
            }
        }
        default:
            Rcpp::stop("Unsupported forest type: expected either a `SampleForest`"
                       " or a `PhylogeneticForest` object.");
    }
}

SEXP get_genome_tour(const SEXP& forest, const bool only_leaves)
{
    using namespace Rcpp;

    switch (TYPEOF(forest)) {
        case S4SXP:
        {
            S4 s4obj(forest);

            if (s4obj.is("Rcpp_PhylogeneticForest")) {
                Environment env(s4obj);

                XPtr<PhylogeneticForest> forest_ptr(env.get(".pointer"));


                LabelTour<PhylogeneticForest> label_tour{*forest_ptr, only_leaves};

                return wrap(label_tour);
            }
        }
        default:
            Rcpp::stop("The first parameter must be a `PhylogeneticForest`.");
    }
}
