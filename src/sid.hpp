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

#ifndef __PROCESS_SID__
#define __PROCESS_SID__

#include <string>

#include <sid.hpp>
#include <mutation_spec.hpp>

#include <Rcpp.h>

class SIDMut : public RACES::Mutations::MutationSpec<RACES::Mutations::SID>
{

    SIDMut(const RACES::Mutations::ChromosomeId& chromosome_id,
           const RACES::Mutations::ChrPosition& chromosomic_position,
           const RACES::Mutations::AlleleId allele_id,
           const std::string& ref, const std::string& alt,
           const std::string& cause="");

public:
    SIDMut();

    inline std::string get_chromosome() const
    {
        return RACES::Mutations::GenomicPosition::chrtos(chr_id);
    }

    inline const RACES::Mutations::ChrPosition& get_position_in_chromosome() const
    {
        return position;
    }

    inline std::string get_ref() const
    {
        return ref;
    }

    inline std::string get_alt() const
    {
        return alt;
    }

    SEXP get_cause() const;

    Rcpp::List get_dataframe() const;

    void show() const;

    static
    SIDMut build_SNV(const SEXP chromosome_name,
                     const SEXP position_in_chromosome,
                     const SEXP alt_base, const SEXP ref_base,
                     const SEXP allele_id, const SEXP cause);

    static
    SIDMut build_SID(const SEXP chromosome_name,
                     const SEXP position_in_chromosome,
                     const SEXP alt_base, const SEXP ref_base,
                     const SEXP allele_id, const SEXP cause);
};

RCPP_EXPOSED_CLASS(SIDMut);

#endif // __PROCESS_SID__
