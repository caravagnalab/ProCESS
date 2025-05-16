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

#include "sid.hpp"

#include "utility.hpp"

SIDMut::SIDMut(const RACES::Mutations::ChromosomeId& chromosome_id,
               const RACES::Mutations::ChrPosition& chromosomic_position,
               const RACES::Mutations::AlleleId allele_id,
               const std::string& ref, const std::string& alt,
               const std::string& cause):
    RACES::Mutations::MutationSpec<RACES::Mutations::SID>(allele_id, chromosome_id,
                                                          chromosomic_position,
                                                          ref, alt, cause)
{}

std::string alleletostr(const RACES::Mutations::AlleleId& allele_id)
{
    if (allele_id == RANDOM_ALLELE) {
        return "random";
    }

    return std::to_string(allele_id);
}

SEXP wrap_allele(const RACES::Mutations::AlleleId& allele_id)
{
    Rcpp::StringVector allele_v(1);

    allele_v[0] = alleletostr(allele_id);

    return allele_v;
}

SIDMut::SIDMut()
{}

SEXP SIDMut::get_cause() const
{
    Rcpp::StringVector cause_v(1);

    if (cause=="") {
        cause_v[0] = NA_STRING;
    } else {
        cause_v[0] = cause;
    }

    return cause_v;
}

Rcpp::List SIDMut::get_dataframe() const
{
    using namespace Rcpp;
    using namespace RACES::Mutations;

    return DataFrame::create(_["chr"]=get_chromosome(),
                             _["chr_pos"]=position,
                             _["allele"]=wrap_allele(allele_id),
                             _["ref"]=get_ref(),
                             _["alt"]=get_alt(),
                             _["type"]=(is_SBS()?"SNV":"indel"),
                             _["cause"]=get_cause());
}

void SIDMut::show() const
{
    using namespace Rcpp;

    if (is_SBS()) {
        Rcout << "SNV";
    } else {
        Rcout << "indel";
    }

    Rcout << "(chr: "<< get_chromosome()
          << ", chr_pos: " << static_cast<size_t>(position)
          << ", allele: " << alleletostr(allele_id)
          << ", ref: " << (ref.size()==0?"-":ref)
          << ", alt: " << (alt.size()==0?"-":alt);

    if (cause!="") {
        Rcout << ", cause: \"" << cause << "\"";
    }
    Rcout << ")" << std::endl;
}

SIDMut SIDMut::build_SNV(const SEXP chromosome_name,
                         const SEXP position_in_chromosome,
                         const SEXP alt_base, const SEXP ref_base,
                         const SEXP allele_id, const SEXP cause)
{
    SIDMut mutation = build_SID(chromosome_name, position_in_chromosome,
                                ref_base, alt_base, allele_id, cause);

    if (mutation.ref.size() != 1) {
        throw std::domain_error("The reference base must be a single "
                                "nucleotide.");
    }

    if (mutation.alt.size() != 1) {
        throw std::domain_error("The altered base must be a single "
                                "nucleotide.");
    }

    return mutation;
}

SIDMut SIDMut::build_SID(const SEXP chromosome_name,
                         const SEXP position_in_chromosome,
                         const SEXP ref_base, const SEXP alt_base,
                         const SEXP allele_id, const SEXP cause)
{
    using namespace Rcpp;
    using namespace RACES::Mutations;

    auto chr_name = as<std::string>(chromosome_name);
    auto chr_id = GenomicPosition::stochr(chr_name);

    auto pos = Rcpp::as<long int>(position_in_chromosome);
    if (pos < 0) {
        throw std::domain_error("Position in chromosome must be a "
                                "non-negative number");
    }

    auto ref_base_str = Rcpp::as<std::string>(ref_base);
    auto alt_base_str = Rcpp::as<std::string>(alt_base);

    auto cause_str = Rcpp::as<std::string>(cause);
    return SIDMut(chr_id, static_cast<RACES::Mutations::ChrPosition>(pos),
                  get_allele_id(allele_id, "allele"), ref_base_str,
                  alt_base_str, cause_str);
}