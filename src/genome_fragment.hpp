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

#ifndef __PROCESS_GENOME_FRAGMENT__
#define __PROCESS_GENOME_FRAGMENT__

#include <string>

#include <genome_fragment.hpp>

#include <Rcpp.h>

class GenomeFragment : private CLONES::Mutations::GenomeFragment
{
    size_t allele_id;

  public:
    GenomeFragment(const std::string& reference_fragment,
                   const size_t& fragment_offset,
                   const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& germline,
                   const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& somatic,
                   const size_t allele_id,
                   const CLONES::Mutations::GenomicPosition& begin_pos,
                   const size_t& size);

    GenomeFragment(const std::string& reference,
                   const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& germline,
                   const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& somatic,
                   const size_t allele_id,
                   const CLONES::Mutations::GenomicPosition& begin_pos,
                   const size_t& size);

    inline std::string get_CIGAR() const
    {
        return CLONES::Mutations::GenomeFragment::get_CIGAR();
    }

    inline const std::string& get_sequence() const
    {
        return CLONES::Mutations::GenomeFragment::get_sequence();
    }

    Rcpp::List get_covered_reference_region() const;

    inline size_t size() const
    {
        return CLONES::Mutations::GenomeFragment::size();
    }

    Rcpp::DataFrame get_mutations() const;

    void show() const;
};

RCPP_EXPOSED_CLASS(GenomeFragment)

#endif // __PROCESS_GENOME_FRAGMENT__
