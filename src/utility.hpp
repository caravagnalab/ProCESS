/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2024 Alberto Casagrande <alberto.casagrande@uniud.it>
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

#ifndef __PROCESS_UTILITY__
#define __PROCESS_UTILITY__

#include <string>
#include <filesystem>

#include <Rcpp.h>

#include <allele.hpp>
#include <archive.hpp>

std::filesystem::path get_tmp_dir_path(const std::string& base_name="ProCESS");

RACES::Mutations::AlleleId get_allele_id(const SEXP allele_id,
                                         const std::string& parameter_name);

std::string ordinal_suffix(const size_t& ord);

inline std::string ordtostr(const size_t ord)
{
    return std::to_string(ord) + ordinal_suffix(ord);
}

template<typename RESULT_TYPE = int,
         std::enable_if_t<std::is_integral_v<RESULT_TYPE>, bool> = true>
RESULT_TYPE get_random_seed(const SEXP seed)
{
    switch (TYPEOF(seed)) {
        case INTSXP:
        case REALSXP:
            return Rcpp::as<RESULT_TYPE>(seed);
        case NILSXP:
        {
            GetRNGstate();
            auto r_seed = R::runif(std::numeric_limits<RESULT_TYPE>::min(),
                                   std::numeric_limits<RESULT_TYPE>::max());

            PutRNGstate();
            return static_cast<RESULT_TYPE>(r_seed);
        }
        default:
            break;
    }

    throw std::domain_error("The seed must be either a "
                            "number or NILL.");
}

void raise_error(const RACES::Archive::WrongFileFormatDescr& exception,
                 const std::string& file_description);

void raise_error(const RACES::Archive::WrongFileFormatVersion& exception,
                 const std::string& file_description);

#endif // __PROCESS_UTILITY__
