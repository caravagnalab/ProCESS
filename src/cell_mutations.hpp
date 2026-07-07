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

#ifndef __PROCESS_CELL_MUTATIONS__
#define __PROCESS_CELL_MUTATIONS__

#include <Rcpp.h>

#include <genome_mutations.hpp>

#include "genome_fragment.hpp"

class GenomeMutations
{
    std::filesystem::path reference_path;
    CLONES::Mutations::GenomeMutations germline;
    CLONES::Mutations::GenomeMutations somatic;

    static size_t count_mutations(const CLONES::Mutations::GenomeMutations& mutations);

    static size_t fill_df(Rcpp::DataFrame df, size_t from,
                          const CLONES::Mutations::GenomeMutations& mutations);
public:
    GenomeMutations();

    GenomeMutations(const std::filesystem::path& reference_path,
                    const CLONES::Mutations::GenomeMutations& germline);

    inline Rcpp::DataFrame get_mutations() const
    {
        return get_mutations(true);
    }

    Rcpp::DataFrame get_mutations(const bool with_germline) const;

    bool apply_contained(const CLONES::Mutations::MutationList& mutation_list);

    GenomeFragment get_fragment(const std::string& chromosome_name,
                                const size_t& allele_id,
                                const size_t& from,
                                const size_t& size) const;

    GenomeFragment get_fragment_from_ref(const std::string& reference_fragment,
                                         const size_t& fragment_offset,
                                         const std::string& chromosome_name,
                                         const size_t& allele_id,
                                         const size_t& from, const size_t& size) const;

    Rcpp::DataFrame get_allele_fragments() const;

    void show() const;
};

RCPP_EXPOSED_CLASS(GenomeMutations)

#endif // __PROCESS_CELL_MUTATIONS__
