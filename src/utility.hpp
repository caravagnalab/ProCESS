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

#ifndef __PROCESS_UTILITY__
#define __PROCESS_UTILITY__

#include <filesystem>
#include <string>
#include <sstream>
#include <concepts>

#include <Rcpp.h>

#include <allele.hpp>
#include <archive.hpp>

std::filesystem::path get_tmp_dir_path(const std::string &base_name = "ProCESS");

CLONES::Mutations::AlleleId get_allele_id(const SEXP allele_id,
                                          const std::string &parameter_name);

class FromSEXP
{
    static bool is_vowel(const char& c);

    template <typename T> struct r_type_traits;

    template <typename T>
     requires std::integral<T> && (!std::same_as<T, bool>)
    struct r_type_traits<T>
    {
        static const int sexp_type = INTSXP;
    };

    template <typename T>
     requires std::floating_point<T>
    struct r_type_traits<T>
    {
        static const int sexp_type = REALSXP;
    };

    template <typename T>
     requires std::is_same_v<T, bool>
    struct r_type_traits<T>
    {
        static const int sexp_type = LGLSXP;
    };


    template <typename T>
     requires std::is_same_v<T, std::string>
    struct r_type_traits<T>
    {
        static const int sexp_type = STRSXP;
    };

    template <typename T>
     requires std::is_same_v<T, Rcpp::List>
    struct r_type_traits<T>
    {
        static constexpr int sexp_type = VECSXP;
    };

public:

    template<typename T>
    static T get(const SEXP& obj, const char* what, const char* expected_type_descr)
    {
        bool error{false};
        if constexpr (std::integral<T> && (!std::same_as<T, bool>)) {
            if (TYPEOF(obj) == REALSXP) {
                double value = Rcpp::as<double>(obj);

                if (std::trunc(value) != value) {
                    std::ostringstream oss;

                    oss << "The " << what << " must be "
                        << (is_vowel(expected_type_descr[0]) ? "an " : "a ")
                        << expected_type_descr << ". A decimal value has"
                        << " been provided (" << value << ").";
                    Rcpp::stop(oss.str());
                }
            }

            error = (TYPEOF(obj) != REALSXP && TYPEOF(obj) != INTSXP);
        } else {
            error = (TYPEOF(obj) != r_type_traits<T>::sexp_type);
        }

        if (error) {
            std::ostringstream oss;
            std::string actual_typename{Rf_type2char(TYPEOF(obj))};

            oss << "The " << what << " must be "
                << (is_vowel(expected_type_descr[0]) ? "an " : "a ")
                << expected_type_descr << ". "
                << (is_vowel(actual_typename[0]) ? "An " : "A ")
                << actual_typename << " has been provided.";

            Rcpp::stop(oss.str());
        }

        return Rcpp::as<T>(obj);
    }
};

template<typename TYPE>
  requires(std::is_integral_v<TYPE>)
std::string ordinal_suffix(const TYPE number) {

    if (number % 100 >= 11 && number % 100 <= 13) {
        return "th";
    }

    // Check the last digit
    switch (number % 10) {
        case 1:  return "st";
        case 2:  return "nd";
        case 3:  return "rd";
        default: return "th";
    }
}

inline std::string ordtostr(const size_t ord)
{
    return std::to_string(ord) + ordinal_suffix(ord);
}

template <typename RESULT_TYPE = int,
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

    Rcpp::stop("The seed must be either a "
               "number or NILL.");
}

void raise_error(const CLONES::Archive::WrongFileFormatDescr &exception,
                 const std::string &file_description);

void raise_error(const CLONES::Archive::WrongFileFormatVersion &exception,
                 const std::string &file_description);

std::string get_demangled_type_name(const std::type_info& type_info);

#endif // __PROCESS_UTILITY__
