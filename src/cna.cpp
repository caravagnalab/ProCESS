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

#include "cna.hpp"

#include "utility.hpp"

// amplification
CNA::CNA(const CLONES::Mutations::GenomicPosition &initial_position,
         const CLONES::Mutations::CNA::Length &length,
         const CLONES::Mutations::AlleleId &allele,
         const CLONES::Mutations::AlleleId &src_allele)
    : CLONES::Mutations::CNA(initial_position, length,
                            CLONES::Mutations::CNA::Type::AMPLIFICATION, src_allele,
                            allele)
{}

// deletion
CNA::CNA(const CLONES::Mutations::GenomicPosition &initial_position,
         const CLONES::Mutations::CNA::Length &length,
         const CLONES::Mutations::AlleleId &allele)
    : CLONES::Mutations::CNA(initial_position, length,
                             CLONES::Mutations::CNA::Type::DELETION, RANDOM_ALLELE,
                             allele)
{}

CNA::CNA() {}

SEXP wrap_allele_id(const CLONES::Mutations::AlleleId &allele_id)
{
    if (allele_id == RANDOM_ALLELE) {
        return Rcpp::wrap(NA_INTEGER);
    }
    return Rcpp::wrap(allele_id);
}

SEXP CNA::get_src_allele() const
{
    if (type == CLONES::Mutations::CNA::Type::AMPLIFICATION) {
        return wrap_allele_id(source);
    }
    return Rcpp::wrap(NA_INTEGER);
}

SEXP CNA::get_allele() const { return wrap_allele_id(dest); }

Rcpp::List CNA::get_dataframe() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    return DataFrame::create(_["chr"] = get_chromosome(),
                             _["from"] = get_position_in_chromosome(),
                             _["length"] = get_length(), _["allele"] = get_allele(),
                             _["src_allele"] = get_src_allele(), _["type"] = get_type());
}

void CNA::show() const
{
    using namespace Rcpp;

    Rcout << "CNA(type: \"" << get_type() << "\", chr: \"" << get_chromosome()
          << "\", from: " << get_position_in_chromosome()
          << ", length: " << get_length();

    if (dest != RANDOM_ALLELE) {
        Rcout << ", allele: " << dest;
    }

    if (get_type() == "A") {
        if (dest != RANDOM_ALLELE) {
            Rcout << ", src_allele: " << source;
        }
    }

    Rcout << ")" << std::endl;
}

CNA CNA::build_CNA(const std::string type, const SEXP chromosome, const SEXP pos_in_chr,
                   const SEXP length, const SEXP allele, const SEXP src_allele)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    auto chr_name = as<std::string>(chromosome);
    auto chr_id = GenomicPosition::stochr(chr_name);

    auto pos = Rcpp::as<long int>(pos_in_chr);
    if (pos < 0) {
        Rcpp::stop("Positions in chromosome must be non-negative numbers.");
    }

    GenomicPosition gen_pos(chr_id, pos);

    auto len = Rcpp::as<long int>(length);
    if (len < 0) {
        Rcpp::stop("Region lengths must be non-negative numbers.");
    }

    AlleleId allele_id = get_allele_id(allele, "allele");

    if (type == "D") {
        return CNA(gen_pos, len, allele_id);
    }

    if (type == "A") {
        AlleleId src_allele_id = get_allele_id(src_allele, "src_allele");

        return CNA(gen_pos, len, allele_id, src_allele_id);
    }

    Rcpp::stop("Unknown CNA type \"" + type +
               "\". Supported types are \"A\" and \"D\" for " +
               "amplification and deletion, respectively");
}
