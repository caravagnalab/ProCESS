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

#include "forest.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"


SEXP ForestCore::load_forest(const std::string &filename, const bool quiet)
{
    if (!std::filesystem::exists(filename)) {
        Rcpp::stop("The file \"" + filename + "\" does not exist.");
    }

    if (!std::filesystem::is_regular_file(filename)) {
        Rcpp::stop("The file \"" + filename + "\" is not a regular file.");
    }

    std::string file_format_description;
    uint8_t file_format_version;
    { // the archive is only used to read the header file
        CLONES::Archive::Binary::In in_archive(filename);

        CLONES::Archive::Binary::In::read_header(in_archive, file_format_description,
                                                 file_format_version);
    }

    if (file_format_description.starts_with(PhylogeneticForest::file_format_header)) {
        return Rcpp::wrap(PhylogeneticForest::load(filename, quiet));
    }

    if (file_format_description.starts_with(SampleForest::file_format_header)) {
        return Rcpp::wrap(SampleForest::load(filename, quiet));
    }

    Rcpp::stop("\"" + filename +"\" is not in a forest file produced by ProCESS.");
}
