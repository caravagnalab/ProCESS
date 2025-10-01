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

#include "sample_forest.hpp"

#include "simulation.hpp"
#include "utility.hpp"

SampleForest::SampleForest() : RACES::Mutants::DescendantForest() {}

SampleForest::SampleForest(const RACES::Mutants::Evolutions::Simulation &simulation)
    : RACES::Mutants::DescendantForest(simulation)
{}

Rcpp::List SampleForest::get_samples_info() const
{
    return TissueSimulation::get_samples_info(get_samples());
}

SampleForest
SampleForest::get_subforest_for(const std::vector<std::string> &sample_names) const
{
    SampleForest forest;

    static_cast<RACES::Mutants::DescendantForest &>(forest) =
        RACES::Mutants::DescendantForest::get_subforest_for(sample_names);

    return forest;
}

void SampleForest::save(const std::string &filename) const
{
    RACES::Archive::Binary::Out out_archive(filename);

    RACES::Mutants::DescendantForest::save(out_archive);
}

SampleForest SampleForest::load(const std::string &filename)
{
    SampleForest forest;

    if (!std::filesystem::exists(filename)) {
        Rcpp::stop("The file \"" + filename + "\" does not exist.");
    }

    if (!std::filesystem::is_regular_file(filename)) {
        Rcpp::stop("The file \"" + filename + "\" is not a regular file.");
    }

    RACES::Archive::Binary::In in_archive(filename);

    try {
        static_cast<RACES::Mutants::DescendantForest &>(forest) =
            RACES::Mutants::DescendantForest::load(in_archive);
    } catch (RACES::Archive::WrongFileFormatDescr &ex) {
        raise_error(ex, "sample forest");
    } catch (RACES::Archive::WrongFileFormatVersion &ex) {
        raise_error(ex, "sample forest");
    }

    return forest;
}

void SampleForest::show() const
{
    using namespace Rcpp;

    size_t num_of_leaves{0};
    for (const auto &sample : get_samples()) {
        num_of_leaves += sample.get_cell_ids().size();
    }

    Rcout << "SampleForest" << std::endl
          << "  # of trees: " << get_roots().size() << std::endl
          << "  # of nodes: " << num_of_nodes() << std::endl
          << "  # of leaves: " << num_of_leaves << std::endl
          << "  samples: {";

    std::string sep = "";
    for (const auto &sample : get_samples()) {
        Rcout << sep << "\"" << sample.get_name() << "\"";
        sep = ", ";
    }

    Rcout << "}" << std::endl << std::endl;
}
