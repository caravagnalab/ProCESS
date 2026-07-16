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

#include "forest.hpp"
#include "sample_forest.hpp"
#include "phylogenetic_forest.hpp"


SEXP ForestCore::load_forest(const std::string &filename, const bool quiet)
{
    if (!std::filesystem::exists(filename)) {
        Rcpp::stop("The file \"" + filename + "\" does not exist.");
    }

    if (!std::filesystem::is_regular_file(filename)) {
        Rcpp::stop("The file \"" + filename + "\" is not a regular file.");
    }

    std::string file_format_description;
    uint8_t file_format_version;
    { // the archive is only used to read the header file
        CLONES::Archive::Binary::In in_archive(filename);

        CLONES::Archive::Binary::In::read_header(in_archive, file_format_description,
                                                 file_format_version);
    }

    if (file_format_description.starts_with(PhylogeneticForest::file_format_header)) {
        return Rcpp::wrap(PhylogeneticForest::load(filename, quiet));
    }

    if (file_format_description.starts_with(SampleForest::file_format_header)) {
        return Rcpp::wrap(SampleForest::load(filename, quiet));
    }

    Rcpp::stop("\"" + filename +"\" is not in a forest file produced by ProCESS.");
}

inline std::string find_code(const CLONES::Mutations::MutationSpec<CLONES::Mutations::SID> &sid,
                             const std::map<CLONES::Mutations::SID, std::string> &driver_codes)
{
    const auto found = driver_codes.find(static_cast<CLONES::Mutations::SID>(sid));
    if (found == driver_codes.end()) {
        return "";
    }

    return (found->second);
}

bool select_column(Rcpp::CharacterVector& column, Rcpp::DataFrame &df, const char* name)
{
    if (df.containsElementNamed(name)) {
        column = df[name];

        return true;
    }

    return false;
}

inline void fill_SID_row(Rcpp::DataFrame& df, 
                         const CLONES::Mutations::MutationSpec<CLONES::Mutations::SID> &sid,
                         const std::map<CLONES::Mutations::SID, std::string> &driver_codes,
                         const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    CharacterVector chrs = df["chr"], types = df["type"], CNA_types = df["CNA_type"],
                    refs = df["ref"], alts = df["alt"], codes = df["code"],
                    natures, causes;

    const bool has_natures = select_column(natures, df, "nature");
    const bool has_causes = select_column(causes, df, "cause");

    IntegerVector starts = df["start"], ends = df["end"], alleles = df["allele"],
                  src_alleles = df["src.allele"];
    
    types[i] = "SID";
    CNA_types[i] = NA_STRING;
    refs[i] = sid.ref;
    alts[i] = sid.alt;
    chrs[i] = GenomicPosition::chrtos(sid.chr_id);
    starts[i] = sid.position;
    ends[i] = sid.position + sid.ref.size() - 1;
    alleles[i] = sid.allele_id;
    src_alleles[i] = NA_INTEGER;

    if (has_natures) {
        natures[i] = sid.get_nature_description();
    }

    if (has_causes) {
        causes[i] = sid.cause;
    }

    const auto code = find_code(sid, driver_codes);

    if (code == "") {
        codes[i] = NA_STRING;
    } else {
        codes[i] = code;
    }
}

inline void fill_CNA_row(Rcpp::DataFrame &df, const CLONES::Mutations::CNA &cna, const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    CharacterVector chrs = df["chr"], types = df["type"], CNA_types = df["CNA_type"],
                    refs = df["ref"], alts = df["alt"], codes = df["code"],
                    natures, causes;

    const bool has_natures = select_column(natures, df, "nature");
    const bool has_causes = select_column(causes, df, "cause");

    IntegerVector starts = df["start"], ends = df["end"], alleles = df["allele"],
                  src_alleles = df["src.allele"];

    types[i] = "CNA";
    refs[i] = NA_STRING;
    alts[i] = NA_STRING;
    chrs[i] = GenomicPosition::chrtos(cna.chr_id);
    starts[i] = cna.position;
    ends[i] = cna.position + cna.length - 1;
    alleles[i] = (cna.dest == RANDOM_ALLELE ? NA_INTEGER : cna.dest);

    if (has_natures) {
        natures[i] = cna.get_nature_description();
    }

    if (has_causes) {
        causes[i] = cna.cause;
    }

    if (cna.type == CLONES::Mutations::CNA::Type::AMPLIFICATION) {
        src_alleles[i] = cna.source;
        CNA_types[i] = "A";
    } else {
        src_alleles[i] = NA_INTEGER;
        CNA_types[i] = "D";
    }
    codes[i] = NA_STRING;
}

inline void fill_WGD_row(Rcpp::DataFrame &df, const std::string& mutant_name, const size_t &i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    CharacterVector chrs = df["chr"], types = df["type"], CNA_types = df["CNA_type"],
                    refs = df["ref"], alts = df["alt"], codes = df["code"],
                    natures, causes;

    const bool has_natures = select_column(natures, df, "nature");
    const bool has_causes = select_column(causes, df, "cause");

    IntegerVector starts = df["start"], ends = df["end"], alleles = df["allele"],
                  src_alleles = df["src.allele"];

    types[i] = "WGD";
    CNA_types[i] = NA_STRING;
    refs[i] = NA_STRING;
    alts[i] = NA_STRING;
    chrs[i] = NA_INTEGER;
    starts[i] = NA_INTEGER;
    ends[i] = NA_INTEGER;
    alleles[i] = NA_INTEGER;
    src_alleles[i] = NA_INTEGER;

    if (has_natures) {
        natures[i] = CLONES::Mutations::Mutation::get_nature_description(CLONES::Mutations::Mutation::Nature::DRIVER);
    }

    if (has_causes) {
        causes[i] = mutant_name;
    }

    codes[i] = NA_STRING;
}

void ForestCore::fill_mutation_list(Rcpp::DataFrame &df, const CLONES::Mutations::MutationList& mutations,
                                    const std::string& mutant_name,
                                    const std::map<CLONES::Mutations::SID, std::string>& driver_codes,
                                    size_t &i, size_t &mut_i)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    IntegerVector application_order = df["order"];
    for (auto dm_it = mutations.begin(); dm_it != mutations.end(); ++dm_it, ++mut_i, ++i) {
        application_order[i] = mut_i;
        switch (dm_it.get_type()) {
        case MutationList::MutationType::SID_TURN:
            fill_SID_row(df, dm_it.get_last_SID(), driver_codes, i);
            break;
        case MutationList::MutationType::CNA_TURN:
            fill_CNA_row(df, dm_it.get_last_CNA(), i);
            break;
        case MutationList::MutationType::WGD_TURN:
            fill_WGD_row(df, mutant_name, i);
            break;
        default:
            Rcpp::stop("Unsupported mutation type " +
                       std::to_string(dm_it.get_type()));
        }
    }
}