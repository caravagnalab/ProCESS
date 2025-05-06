/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2025 Alberto Casagrande <alberto.casagrande@uniud.it>
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

#ifndef __PROCESS_SEQ_SIMULATION__
#define __PROCESS_SEQ_SIMULATION__

#include <Rcpp.h>

#include <read_simulator.hpp>

#include "sequencers.hpp"
#include "phylogenetic_forest.hpp"

Rcpp::List  simulate_seq(const PhylogeneticForest& forest, SEXP& sequencer,
                         SEXP& reference_genome,
                         SEXP& chromosome_ids, const double& coverage,
                         const int& read_size, const int& insert_size_mean,
                         const int& insert_size_stddev,
                         const std::string& output_dir, const bool& write_SAM,
                         const bool& update_SAM_dir,
                         const SEXP& FACS_labelling_function,
                         const double& purity, const bool& with_normal_sample,
                         const bool& preneoplastic_in_normal,
                         const std::string& filename_prefix,
                         const std::string& template_name_prefix,
                         const bool& include_non_sequenced_mutations,
                         const SEXP& seed);

Rcpp::List  simulate_normal_seq(const PhylogeneticForest& forest, SEXP& sequencer,
                                SEXP& reference_genome,
                                SEXP& chromosome_ids, const double& coverage,
                                const int& read_size, const int& insert_size_mean,
                                const int& insert_size_stddev,
                                const std::string& output_dir, const bool& write_SAM,
                                const bool& update_SAM_dir,
                                const bool& with_preneoplastic,
                                const std::string& filename_prefix,
                                const std::string& template_name_prefix,
                                const bool& include_non_sequenced_mutations,
                                const SEXP& seed);

#endif // __PROCESS_SEQ_SIMULATION__
