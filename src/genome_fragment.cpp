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

#include "genome_fragment.hpp"

GenomeFragment::GenomeFragment(const std::string& reference_fragment,
                               const size_t& fragment_offset,
                               const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& germline,
                               const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& somatic,
                               const size_t allele_id,
                               const CLONES::Mutations::GenomicPosition& begin_pos,
                               const size_t& size):
    CLONES::Mutations::GenomeFragment{reference_fragment, fragment_offset, germline,
                                       somatic, begin_pos, size},
    allele_id{allele_id}
{}

GenomeFragment::GenomeFragment(const std::string& reference,
                               const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& germline,
                               const std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>>& somatic,
                               const size_t allele_id,
                               const CLONES::Mutations::GenomicPosition& begin_pos,
                               const size_t& size):
    GenomeFragment{reference, 0, germline, somatic, allele_id, begin_pos, size}
{}

Rcpp::List GenomeFragment::get_covered_reference_region() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    auto g_region = CLONES::Mutations::GenomeFragment::get_covered_reference_region();

    return List::create(_["chr"] = GenomicPosition::chrtos(g_region.get_chromosome_id()),
                        _["allele"] = allele_id,
                        _["from"] = g_region.begin(),
                        _["size"] = g_region.size());
}

Rcpp::DataFrame GenomeFragment::get_mutations() const
{
    using namespace Rcpp;

    const auto& mutations = CLONES::Mutations::GenomeFragment::get_mutations();

    const size_t nrows{mutations.size()};

    IntegerVector chr_pos(nrows), alleles(nrows);
    CharacterVector chr_names(nrows), refs(nrows),
        alts(nrows), causes(nrows), classes(nrows);

    std::string chr_name;
    if (nrows > 0) {
        chr_name = CLONES::Mutations::GenomicPosition::chrtos(mutations.begin()->chr_id);
    }

    size_t i{0};
    for (const auto& sid : mutations) {
        chr_names[i] = chr_name;
        alleles[i] = allele_id;
        chr_pos[i] = sid.position;
        refs[i] = sid.ref;
        alts[i] = sid.alt;
        causes[i] = sid.cause;
        classes[i] = CLONES::Mutations::Mutation::get_nature_description(sid.nature);

        ++i;
    }

    return DataFrame::create(_["chr"] = chr_names, _["allele"] = alleles,
                             _["from"] = chr_pos,
                             _["ref"] = refs, _["alt"] = alts,
                             _["causes"] = causes, _["classes"] = classes);
}

void GenomeFragment::show() const
{
    using namespace Rcpp;

    const auto g_pos{CLONES::Mutations::GenomeFragment::get_genomic_position()};

    const auto from = static_cast<size_t>(g_pos.position);

    Rcout << "chr" << CLONES::Mutations::GenomicPosition::chrtos(g_pos.chr_id)
          << "(" << static_cast<size_t>(allele_id) << ")["
          << from << "-"
          << (from + CLONES::Mutations::GenomeFragment::size()) << "]";
}