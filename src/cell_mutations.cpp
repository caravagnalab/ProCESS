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

#include <sstream>

#include <Rcpp.h>

#include <fasta_chr_reader.hpp>

#include "cell_mutations.hpp"


GenomeMutations::GenomeMutations()
{}

GenomeMutations::GenomeMutations(const std::filesystem::path& reference_path,
                                 const CLONES::Mutations::GenomeMutations& germline):
    reference_path{reference_path}, germline{germline}, somatic{germline.copy_structure()}
{}

size_t GenomeMutations::count_mutations(const CLONES::Mutations::GenomeMutations& mutations)
{
    size_t num_of_mutations{0};
    for (const auto& [chr_id, chr_mutations] : mutations.get_chromosomes()) {
        for (const auto& [allele_id, allele_mutations] : chr_mutations.get_alleles()) {
            for (const auto& [pos, fragment_mutations]: allele_mutations.get_fragments()) {
                num_of_mutations += fragment_mutations.get_mutations().size();
            }
        }
    }

    return num_of_mutations;
}

size_t GenomeMutations::fill_df(Rcpp::DataFrame df, size_t from,
                                const CLONES::Mutations::GenomeMutations& mutations)
{
    using namespace Rcpp;

    IntegerVector chr_pos = df["from"], alleles = df["allele"];
    CharacterVector chr_names = df["chr"], refs = df["ref"],
        alts = df["alt"], causes = df["causes"], classes = df["classes"];

    for (const auto& [chr_id, chr_mutations] : mutations.get_chromosomes()) {
        const std::string chr_name = CLONES::Mutations::GenomicPosition::chrtos(chr_id);
        for (const auto& [allele_id, allele_mutations] : chr_mutations.get_alleles()) {
            for (const auto& [f_pos, fragment_mutations]: allele_mutations.get_fragments()) {
                for (const auto& [pos, mutation_ptr]: fragment_mutations.get_mutations()) {
                    chr_names[from] = chr_name;
                    chr_pos[from] = pos.position;

                    refs[from] = mutation_ptr->ref;
                    alts[from] = mutation_ptr->alt;
                    causes[from] = mutation_ptr->cause;
                    classes[from] = CLONES::Mutations::Mutation::get_nature_description(mutation_ptr->nature);
                    alleles[from] = allele_id;

                    ++from;
                }
            }
        }
    }

    return from;
}

Rcpp::DataFrame GenomeMutations::get_mutations(const bool with_germline) const
{
    size_t nrows{count_mutations(somatic)};

    if (with_germline) {
        nrows += count_mutations(germline);
    }

    using namespace Rcpp;

    IntegerVector chr_pos(nrows), alleles(nrows);
    CharacterVector chr_names(nrows), refs(nrows),
        alts(nrows), causes(nrows), classes(nrows);

    DataFrame df = DataFrame::create(_["chr"] = chr_names, _["allele"] = alleles,
                                     _["from"] = chr_pos,
                                     _["ref"] = refs, _["alt"] = alts,
                                     _["causes"] = causes, _["classes"] = classes);

    const size_t from = fill_df(df, 0, somatic);

    if (with_germline) {
        fill_df(df, from, germline);
    }

    return df;
}

bool GenomeMutations::apply_contained(const CLONES::Mutations::MutationList& mutation_list)
{
    return somatic.apply_contained(mutation_list);
}

GenomeFragment GenomeMutations::get_fragment(const std::string& chromosome_name,
                                             const size_t& allele_id,
                                             const size_t& from,
                                             const size_t& size) const
{
    using namespace CLONES::Mutations;
    using namespace CLONES::IO;
    using namespace CLONES::IO::FASTA;

    IndexedReader<ChromosomeData<Sequence>> fasta_reader(reference_path);

    const size_t delta{100};

    size_t offset{(from<delta?0: from-delta)};

    std::string reference_fragment;
    fasta_reader.read(reference_fragment, chromosome_name, offset, size+delta+from-offset);

    return get_fragment_from_ref(reference_fragment, offset,
                                 chromosome_name, allele_id, from, size);
}

std::pair<std::map<CLONES::Mutations::GenomicPosition, std::shared_ptr<CLONES::Mutations::SID>> const *,
          CLONES::Mutations::AlleleId>
get_mutations_in_fragment(const CLONES::Mutations::GenomeMutations& mutations,
                          const CLONES::Mutations::ChromosomeId& chr_id,
                          const size_t& allele_id, const size_t& from, const size_t& size)
{
    using namespace CLONES::Mutations;

    const auto& chr_mutations = mutations.get_chromosome(chr_id);

    const auto& allele_mutations = chr_mutations.get_allele(static_cast<AlleleId>(allele_id));

    AlleleId src_allele_id = allele_mutations.get_history().front();

    GenomicPosition begin_pos{chr_id, static_cast<ChrPosition>(from)};
    GenomicRegion fragment_region{begin_pos, static_cast<GenomicRegion::Length>(size)};
    if (!allele_mutations.contains(fragment_region)) {

        std::ostringstream oss;

        oss << "The allele " << std::to_string(allele_id) << " of chromosome "
            << CLONES::Mutations::GenomicPosition::chrtos(chr_id)
            << " does not contain the region ["
            << fragment_region.begin() << "," << fragment_region.end() << "].";

        Rcpp::stop(oss.str());
    }

    auto fragment_it = allele_mutations.get_fragments().upper_bound(begin_pos);

    if (fragment_it != allele_mutations.get_fragments().begin()) {
        --fragment_it;
    }

    return {&(fragment_it->second.get_mutations()), src_allele_id};
}

GenomeFragment GenomeMutations::get_fragment_from_ref(const std::string& reference_fragment,
                                                      const size_t& fragment_offset,
                                                      const std::string& chromosome_name,
                                                      const size_t& allele_id,
                                                      const size_t& from, const size_t& size) const
{
    try {
        using ChromosomeId = CLONES::Mutations::ChromosomeId;
        using GenomicPosition = CLONES::Mutations::GenomicPosition;

        const ChromosomeId chr_id = GenomicPosition::stochr(chromosome_name);

        auto [somatic_muts, src_allele_id] = get_mutations_in_fragment(somatic, chr_id, allele_id,
                                                                       from, size);

        auto [germline_muts, germline_allele_id] = get_mutations_in_fragment(germline, chr_id,
                                                                             src_allele_id,
                                                                             from, size);

        GenomicPosition begin_pos{chr_id, static_cast<CLONES::Mutations::ChrPosition>(from)};

        return GenomeFragment{reference_fragment, fragment_offset, *germline_muts,
                              *somatic_muts, allele_id, begin_pos, size};
    } catch (const std::exception &ex) {
        Rcpp::stop(ex.what());
    }
}

Rcpp::DataFrame GenomeMutations::get_allele_fragments() const
{
    size_t num_of_fragments{0};
    for (const auto& [chr_id, chr_mutation] : somatic.get_chromosomes()) {
        for (const auto& [allele_id, allele_mutation] : chr_mutation.get_alleles()) {
            num_of_fragments += allele_mutation.get_fragments().size();
        }
    }

    using namespace Rcpp;

    IntegerVector chr_pos(num_of_fragments), alleles(num_of_fragments),
                  sizes(num_of_fragments), allele_srcs(num_of_fragments);
    CharacterVector chr_names(num_of_fragments);

    size_t i{0};
    for (const auto& [chr_id, chr_mutation] : somatic.get_chromosomes()) {
        const std::string chr_name = CLONES::Mutations::GenomicPosition::chrtos(chr_id);
        for (const auto& [allele_id, allele_mutation] : chr_mutation.get_alleles()) {
            for (const auto& [g_pos, allele_fragment] : allele_mutation.get_fragments()) {
                chr_names[i] = chr_name;
                alleles[i] = allele_id;
                allele_srcs[i] = allele_mutation.get_history().back();
                chr_pos[i] = allele_fragment.begin();
                sizes[i] = allele_fragment.size();
                ++i;
            }
        }
    }

    return DataFrame::create(_["chr"] = chr_names, _["allele"] = alleles,
                             _["src_allele"] = allele_srcs,
                             _["from"] = chr_pos, _["size"] = sizes);
}

void GenomeMutations::show() const
{
    using namespace Rcpp;

    Rcout << "GenomeMutations: " << somatic.get_chromosomes().size()
          << " chrs ";

    size_t num_alleles{0};
    for (const auto& [chr_id, chr_mutation] : somatic.get_chromosomes()) {
        num_alleles += chr_mutation.get_alleles().size();
    }

    Rcout << num_alleles << " alleles";
}