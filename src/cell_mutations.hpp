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
#include <union_map_proxy.hpp>

#include "genome_fragment.hpp"

class GenomeMutations
{
public:

    GenomeMutations();

    GenomeMutations(const std::filesystem::path& reference_path,
                    const CLONES::Mutations::GenomeMutations& germline);

    GenomeMutations(const std::filesystem::path& reference_path,
                    const CLONES::Mutations::GenomeMutations& germline,
                    const CLONES::Mutations::GenomeMutations& somatic);

    Rcpp::List region_aligning_on_reference(const std::string& chromosome_name,
                                            const size_t allele_id,
                                            const size_t from,
                                            size_t size) const;

    inline Rcpp::DataFrame get_mutations() const
    {
        return get_mutations(true);
    }

    bool apply_contained(const CLONES::Mutations::MutationList& mutation_list);

    std::list<CLONES::Mutations::AlleleId>
    get_alleles_covering_ref_region(const std::string& chromosome_name,
                                    const size_t& from,
                                    const size_t& size) const;

    GenomeFragment get_fragment(const std::string& chromosome_name,
                                const size_t& allele_id,
                                const size_t& from,
                                const size_t& size) const;

    GenomeFragment get_fragment(const std::string& chromosome_name,
                                const size_t& allele_id,
                                const size_t& from, const size_t& size,
                                const std::string& reference_fragment,
                                const size_t& fragment_offset) const;

    Rcpp::DataFrame get_allele_fragments() const;

    size_t num_of_mutations(const bool with_germline = true) const;

    size_t num_of_CNAs() const;

    size_t fill_mutation_df(Rcpp::DataFrame &df, size_t index = 0, const bool with_germline = true) const;

    inline size_t fill_CNA_df(Rcpp::DataFrame &df, size_t index = 0) const
    {
        // all CNAs appear as somatic CNAs
        return fill_CNA_df(df, somatic, index);
    }

    Rcpp::DataFrame get_mutations(const bool with_germline = true) const;

    Rcpp::DataFrame get_CNAs() const;

    void show() const;

    static Rcpp::DataFrame create_mutation_df(const size_t nrows);

    static Rcpp::DataFrame create_CNA_df(const size_t nrows);

    static size_t fill_mutation_df(Rcpp::DataFrame &df, const CLONES::Mutations::GenomeMutations &mutations,
                                   size_t index = 0);

    static size_t fill_CNA_df(Rcpp::DataFrame &df, const CLONES::Mutations::GenomeMutations &mutations,
                              size_t index = 0);
private:
    using mutation_map = CLONES::union_map_proxy<CLONES::Mutations::GenomicPosition,
                                                 std::shared_ptr<CLONES::Mutations::SID>>;

    mutation_map get_fragment_mutations(const CLONES::Mutations::ChromosomeId& chr_id,
                                        const CLONES::Mutations::AlleleId& allele_id,
                                        const size_t& from) const;

    CLONES::Mutations::GenomicRegion
    static ref_aligning_on_region(const CLONES::Mutations::ChromosomeId& chr_id,
                                  const mutation_map& f_mutations,
                                  const size_t from, const size_t size);

    GenomeFragment get_fragment_from_ref(const std::string& reference_fragment,
                                         const size_t& fragment_offset,
                                         const CLONES::Mutations::ChromosomeId chr_id,
                                         const CLONES::Mutations::AlleleId allele_id,
                                         const size_t from, const size_t size) const;

    std::filesystem::path reference_path;
    CLONES::Mutations::GenomeMutations germline;
    CLONES::Mutations::GenomeMutations somatic;
};

RCPP_EXPOSED_CLASS(GenomeMutations)

#endif // __PROCESS_CELL_MUTATIONS__
