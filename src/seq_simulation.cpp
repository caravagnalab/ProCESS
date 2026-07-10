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

#include <utils.hpp>

#include "seq_simulation.hpp"
#include "sequencers.hpp"

#include "utility.hpp"

std::string join(const std::set<std::string> &S, const char &sep = ';')
{
    std::ostringstream oss;

    if (S.size() > 0) {
        auto S_it = S.begin();
        oss << *S_it;
        while (++S_it != S.end()) {
            oss << sep << *S_it;
        }
    }

    return oss.str();
}

std::set<std::string>
get_descriptions(const std::set<CLONES::Mutations::Mutation::Nature> &nature_set)
{
    std::set<std::string> nature_strings;

    for (const auto &nature : nature_set) {
        nature_strings.insert(CLONES::Mutations::Mutation::get_nature_description(nature));
    }

    return nature_strings;
}

void add_SNV_data(
    Rcpp::DataFrame &df,
    const std::map<CLONES::Mutations::SID,
                   CLONES::Mutations::SequencingSimulations::SIDData> &mutations)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    size_t num_of_mutations = mutations.size();

    IntegerVector chr_pos(num_of_mutations);
    CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
        alt(num_of_mutations), causes(num_of_mutations), natures(num_of_mutations);

    size_t index{0};
    for (const auto &[mutation, data] : mutations) {
        chr_names[index] = GenomicPosition::chrtos(mutation.chr_id);
        chr_pos[index] = mutation.position;
        ref[index] = mutation.ref;
        alt[index] = mutation.alt;

        if (data.causes.size() == 0) {
            causes[index] = NA_STRING;
        } else {
            causes[index] = join(data.causes, ';');
        }

        auto descr_set = get_descriptions(data.nature_set);
        natures[index] = join(descr_set, ';');

        ++index;
    }

    df.push_back(chr_names, "chr");
    df.push_back(chr_pos, "from");
    df.push_back(ref, "ref");
    df.push_back(alt, "alt");
    df.push_back(causes, "cause");
    df.push_back(natures, "nature");
}

void add_wide_sample_statistics(
    Rcpp::DataFrame &df,
    const CLONES::Mutations::SequencingSimulations::SampleStatistics &sample_statistics,
    const std::map<CLONES::Mutations::SID,
                   CLONES::Mutations::SequencingSimulations::SIDData> &mutations)
{
    if (df.length() == 0) {
        add_SNV_data(df, mutations);
    }

    size_t num_of_mutations = mutations.size();

    if (num_of_mutations != static_cast<size_t>(df.nrows())) {
        Rcpp::stop("SeqSimResults are not canonical!!!");
    }

    using namespace Rcpp;
    using namespace CLONES::Mutations;

    DoubleVector VAF(num_of_mutations);
    IntegerVector occurrences(num_of_mutations), coverages(num_of_mutations);

    size_t index{0};
    auto coverage_it = sample_statistics.get_coverage().begin();
    std::less<GenomicPosition> come_before;
    for (const auto &[mutation, mutation_data] : mutations) {

        while (come_before(coverage_it->first, mutation)) {
            ++coverage_it;
        }

        coverages[index] = coverage_it->second;

        auto it = sample_statistics.get_data().find(mutation);
        if (it != sample_statistics.get_data().end()) {
            occurrences[index] = it->second.num_of_occurrences;
            VAF[index] =
                static_cast<double>(it->second.num_of_occurrences) / coverage_it->second;
        } else {
            occurrences[index] = 0;
            VAF[index] = 0;
        }

        ++index;
    }

    const auto &sample_name = sample_statistics.get_sample_name();

    df.push_back(occurrences, sample_name + ".NV");
    df.push_back(coverages, sample_name + ".DP");
    df.push_back(VAF, sample_name + ".VAF");
}

std::map<CLONES::Mutations::SID, CLONES::Mutations::SequencingSimulations::SIDData>
get_active_mutations(const CLONES::Mutations::SequencingSimulations::SampleSetStatistics
                         &sample_set_statistics,
                     const bool &include_non_sequenced_mutations)
{
    std::map<CLONES::Mutations::SID, CLONES::Mutations::SequencingSimulations::SIDData>
        active_mutations;

    for (const auto &[sample_name, sample_stats] : sample_set_statistics) {
        for (const auto &[mutation, mutation_data] : sample_stats.get_data()) {
            if (mutation_data.num_of_occurrences > 0 || include_non_sequenced_mutations) {

                auto it = active_mutations.find(mutation);
                if (it == active_mutations.end()) {
                    active_mutations.insert({mutation, mutation_data});
                } else {
                    auto &data = it->second;
                    data.num_of_occurrences += mutation_data.num_of_occurrences;

                    data.causes = get_union(data.causes, mutation_data.causes);
                    data.nature_set =
                        get_union(data.nature_set, mutation_data.nature_set);
                }
            }
        }
    }

    return active_mutations;
}

Rcpp::List
get_wide_dataframe(const CLONES::Mutations::SequencingSimulations::SampleSetStatistics
                         &sample_set_statistics,
                   const bool &include_non_sequenced_mutations)
{
    const auto mutations =
        get_active_mutations(sample_set_statistics, include_non_sequenced_mutations);

    auto df = Rcpp::DataFrame::create();
    for (const auto &[sample_name, sample_stats] : sample_set_statistics) {
        add_wide_sample_statistics(df, sample_stats, mutations);
    }

    return df;
}


Rcpp::List
get_long_dataframe(const CLONES::Mutations::SequencingSimulations::SampleSetStatistics
                    &sample_set_statistics)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    size_t num_of_rows{0};

    for (const auto &[sample_name, sample_stats] : sample_set_statistics) {
        num_of_rows += sample_stats.get_data().size();
    }

    IntegerVector chr_pos(num_of_rows), occurrences(num_of_rows),
                  coverages(num_of_rows);
    DoubleVector VAF(num_of_rows);
    CharacterVector samples(num_of_rows), chr_names(num_of_rows), ref(num_of_rows),
                    alt(num_of_rows), causes(num_of_rows), natures(num_of_rows);

    size_t index{0};
    for (const auto &[sample_name, sample_stats] : sample_set_statistics) {
        for (const auto &[sid, sid_data] : sample_stats.get_data()) {
            samples[index] = sample_name;
            chr_names[index] = GenomicPosition::chrtos(sid.chr_id);
            chr_pos[index] = sid.position;
            ref[index] = sid.ref;
            alt[index] = sid.alt;
            if (sid_data.causes.size() == 0) {
                causes[index] = NA_STRING;
            } else {
                causes[index] = join(sid_data.causes, ';');
            }

            const auto descr_set = get_descriptions(sid_data.nature_set);
            natures[index] = join(descr_set, ';');
            
            occurrences[index] = sid_data.num_of_occurrences;

            const auto coverage = sample_stats.get_coverage(sid);
            coverages[index] = coverage;

            if (coverage != 0) {
                VAF[index] = static_cast<double>(sid_data.num_of_occurrences) / coverage;
            } else {
                VAF[index] = R_NaN;
            }

            ++index;
        }
    }

    return DataFrame::create(_["sample"] = samples,
                             _["chr"] = chr_names, _["from"] = chr_pos,
                             _["ref"] = ref, _["alt"] = alt,
                             _["cause"] = causes,
                             _["nature"] = natures, _["NV"] = occurrences,
                             _["DP"] = coverages, _["VAF"] = VAF);
}

std::filesystem::path get_reference_genome(const PhylogeneticForest &forest,
                                           const SEXP &reference_genome)
{
    switch (TYPEOF(reference_genome)) {
    case NILSXP:
    {
        const auto ref_genome = forest.get_reference_path();
        if (!std::filesystem::exists(ref_genome)) {
            Rcpp::stop("The reference genome file \"" + to_string(ref_genome) +
                       "\" does not exists anymore. " + "Please, re-build the mutation " +
                       "engine or use the parameter " + "\"reference_genome\".");
        }

        return ref_genome;
    }
    case STRSXP:
    {
        const auto ref_genome = Rcpp::as<std::string>(reference_genome);

        if (!std::filesystem::exists(ref_genome)) {
            Rcpp::stop("The reference genome file \"" + to_string(ref_genome) +
                       "\" does not exists.");
        }

        return ref_genome;
    }
    default:
        Rcpp::stop("The parameter \"reference_genome\" must be "
                   "either NULL or a string.");
    }
}

std::set<CLONES::Mutations::ChromosomeId>
get_relevant_chr_set(const PhylogeneticForest& forest, SEXP &chromosome_ids)
{
    using namespace Rcpp;
    using namespace CLONES::Mutations;

    switch (TYPEOF(chromosome_ids)) {
    case NILSXP:
    {
        std::set<CLONES::Mutations::ChromosomeId> chr_ids;

        const auto& germline_mutations = forest.get_germline_mutations();
        const auto chr_view = std::views::keys(germline_mutations.get_chromosomes());

        return std::set<CLONES::Mutations::ChromosomeId>{chr_view.begin(), chr_view.end()};
    }
    case STRSXP:
    {
        std::set<ChromosomeId> chr_ids;

        CharacterVector chr_names{chromosome_ids};

        for (const auto &chr_name : chr_names) {
            chr_ids.insert(GenomicPosition::stochr(as<std::string>(chr_name)));
        }
        return chr_ids;
    }
    case VECSXP:
    {
        std::set<ChromosomeId> chr_ids;
        List chr_names = Rcpp::as<List>(chromosome_ids);

        size_t i{0};
        for (const auto &chr_name : chr_names) {
            ++i;
            if (TYPEOF(chr_name) != STRSXP) {
                Rcpp::stop("Expected a list of string: the " + ordtostr(i) +
                           " element of the list is not " + "a string.");
            }

            Rcpp::CharacterVector name{chr_name};
            if (name.length() > 1) {
                Rcpp::stop("Expected a list of string: the " + ordtostr(i) +
                           " element of the list is not " + "a string.");
            }

            chr_ids.insert(GenomicPosition::stochr(as<std::string>(name)));
        }
        return chr_ids;
    }
    default:
        Rcpp::stop("Unsupported chromosome list type");
    }
}

template <template <class> typename QUALITY_SCORE_MODEL>
inline CLONES::Mutations::SequencingSimulations::SampleSetStatistics
simulate_seq(CLONES::Mutations::SequencingSimulations::ReadSimulator<> &simulator,
             const Rcpp::XPtr<BasicIlluminaSequencer> &R_seq,
             const PhylogeneticForest& forest,
             const std::set<CLONES::Mutations::ChromosomeId> &chromosome_ids,
             const double& coverage, const bool& produce_normal_sample,
             const double& purity, const bool& with_pre_neoplastic,
             const bool& with_germline, const std::string &base_name,
             const bool& missed_SID_statistics, const bool& germinal_statistics,
             std::ostream &progress_bar_stream,
             const int &seed, const bool quiet)
{
    using namespace CLONES::Sequencers;

    Illumina::BasicSequencer<QUALITY_SCORE_MODEL> sequencer(R_seq->get_error_rate(),
                                                            seed);

    return simulator(sequencer, forest, chromosome_ids, coverage, produce_normal_sample,
                     purity, with_pre_neoplastic, with_germline, base_name,
                     missed_SID_statistics, germinal_statistics, progress_bar_stream,
                     quiet);
}

CLONES::Mutations::SequencingSimulations::SampleSetStatistics simulate_seq(
    CLONES::Mutations::SequencingSimulations::ReadSimulator<> &simulator, SEXP &sequencer,
    const PhylogeneticForest& forest,
    const std::set<CLONES::Mutations::ChromosomeId> &chromosome_ids,
    const double& coverage, const bool& produce_normal_sample,
    const double& purity, const bool& with_pre_neoplastic,
    const bool& with_germline, const std::string &base_name,
    const bool& missed_SID_statistics, const bool& germinal_statistics,
    std::ostream &progress_bar_stream,
    const int &seed, const bool quiet)
{
    switch (TYPEOF(sequencer)) {
    case S4SXP:
    {
        Rcpp::S4 s4obj(sequencer);
        if (s4obj.is("Rcpp_BasicIlluminaSequencer")) {
            Rcpp::Environment env(s4obj);

            Rcpp::XPtr<BasicIlluminaSequencer> sequencer_ptr(env.get(".pointer"));

            using namespace CLONES::Sequencers;

            if (sequencer_ptr->producing_random_scores()) {
                return simulate_seq<QualityScoreModel>(
                    simulator, sequencer_ptr, forest, chromosome_ids, coverage,
                    produce_normal_sample, purity, with_pre_neoplastic, with_germline,
                    base_name, missed_SID_statistics, germinal_statistics,
                    progress_bar_stream, seed, quiet);
            } else {
                return simulate_seq<ConstantQualityScoreModel>(
                    simulator, sequencer_ptr, forest, chromosome_ids, coverage,
                    produce_normal_sample, purity, with_pre_neoplastic, with_germline,
                    base_name, missed_SID_statistics, germinal_statistics,
                    progress_bar_stream, seed, quiet);
            }
        }
        if (s4obj.is("Rcpp_ErrorlessIlluminaSequencer")) {
            CLONES::Sequencers::Illumina::ErrorLessSequencer seq;

            return simulator(seq, forest, chromosome_ids, coverage,
                             produce_normal_sample, purity, with_pre_neoplastic,
                             with_germline, base_name, missed_SID_statistics,
                             germinal_statistics, progress_bar_stream, quiet);
        }

        Rcpp::stop("Unsupported sequencer type");
    }
    case NILSXP:
    {
        CLONES::Sequencers::Illumina::ErrorLessSequencer seq;

        return simulator(seq, forest, chromosome_ids, coverage, produce_normal_sample,
                         purity, with_pre_neoplastic, with_germline, base_name,
                         missed_SID_statistics, germinal_statistics,
                         progress_bar_stream, quiet);
    }
    default:
        Rcpp::stop("Unsupported sequencer type");
    }
}

std::binomial_distribution<uint32_t> get_bin_dist(const int &insert_size_mean,
                                                  const int &insert_size_stddev)
{
    double q =
        static_cast<double>(insert_size_stddev * insert_size_stddev) / insert_size_mean;
    double p = 1 - q;
    if (p < 0) {
        Rcpp::stop("The insert size mean (" + std::to_string(insert_size_mean) +
                   ") must" + " be greater than or equal to its variance (" +
                   std::to_string(insert_size_stddev) + "*" +
                   std::to_string(insert_size_stddev) + "=" +
                   std::to_string(insert_size_stddev * insert_size_stddev) + ").\n" +
                   "Set the standard deviation and the variance by using " +
                   "the optional parameter \"insert_size_stddev\".");
    }

    uint32_t t = static_cast<uint32_t>(insert_size_mean / p);

    return std::binomial_distribution<uint32_t>(t, p);
}

template <typename SEQUENCER_CLASS>
Rcpp::List extract_sequencer_data(SEXP &sequencer, const std::string &seq_class_name)
{
    using namespace Rcpp;

    S4 s4obj(sequencer);
    Environment env(s4obj);

    XPtr<SEQUENCER_CLASS> rel_ptr(env.get(".pointer"));

    if constexpr (std::is_base_of_v<SEQUENCER_CLASS, BasicIlluminaSequencer>) {
        return List::create(
            _["name"] = seq_class_name, _["error_rate"] = rel_ptr->get_error_rate(),
            _["random_quality_scores"] = rel_ptr->producing_random_scores());
    }

    if constexpr (std::is_base_of_v<SEQUENCER_CLASS, ErrorlessIlluminaSequencer>) {
        return List::create(_["name"] = seq_class_name,
                            _["error_rate"] = rel_ptr->get_error_rate());
    }

    Rcpp::stop("Unsupported sequencer type");
}

Rcpp::List get_sequencer_data(SEXP &sequencer)
{
    using namespace Rcpp;

    switch (TYPEOF(sequencer)) {
    case NILSXP:
    {
        return sequencer;
    }
    case S4SXP:
    {
        S4 s4obj(sequencer);

        std::string seq_class_name = s4obj.attr("class");
        seq_class_name = seq_class_name.substr(strlen("Rcpp_"));

        if (s4obj.is("Rcpp_BasicIlluminaSequencer")) {
            return extract_sequencer_data<BasicIlluminaSequencer>(sequencer,
                                                                  seq_class_name);
        }
        if (s4obj.is("Rcpp_ErrorlessIlluminaSequencer")) {
            return extract_sequencer_data<ErrorlessIlluminaSequencer>(sequencer,
                                                                      seq_class_name);
        }

        std::ostringstream oss;

        oss << "Unsupported sequencer class \"" << seq_class_name << "\"";

        Rcpp::stop(oss.str());
    }
    default:
    {
        Function stop("stop"), paste0("paste");

        stop(paste0(sequencer, " is not supported as sequencer"));
    }
    }

    Rcpp::stop("Unsupported sequencer type");
}

Rcpp::List simulate_seq(const PhylogeneticForest &forest, SEXP &sequencer,
                        SEXP &reference_genome, SEXP &chromosome_ids,
                        const double &coverage, const int &read_size,
                        const int &insert_size_mean, const int &insert_size_stddev,
                        const std::string &output_dir, const bool &write_SAM,
                        const bool &update_SAM_dir,
                        const double &purity, const bool &with_normal_sample,
                        const bool &pre_neoplastic_in_normal,
                        const std::string &filename_prefix,
                        const std::string &template_name_prefix,
                        const bool &missed_SID_statistics,
                        const bool &germline_statistics,
                        const bool &wide_format, const SEXP &seed,
                        const bool &quiet)
{
    using namespace CLONES::Mutations::SequencingSimulations;


    if (with_normal_sample && !wide_format
            && !(pre_neoplastic_in_normal || germline_statistics)) {
        Rcpp::warning("\"with_normal_sample\"=TRUE, but neither \"pre_neoplastic_in_normal\" "
                      "nor \"germline_statistics\" are TRUE. The output statistics will "
                      "contain any reference to the normal sample.");
    }

    const auto ref_genome = get_reference_genome(forest, reference_genome);

    std::filesystem::path output_path = output_dir;

    bool remove_output_path = false;

    if (!write_SAM) {
        remove_output_path = true;
        output_path = get_tmp_dir_path(output_dir);
    }

    ReadSimulator<>::Mode SAM_mode = ReadSimulator<>::Mode::CREATE;

    if (update_SAM_dir) {
        SAM_mode = ReadSimulator<>::Mode::UPDATE;
    }

    ReadSimulator<> simulator;
    auto c_seed = get_random_seed<int>(seed);
    if (insert_size_mean == 0) {
        simulator = ReadSimulator<>(output_path, ref_genome, read_size, SAM_mode, false,
                                    template_name_prefix, c_seed);
    } else {
        auto insert_size_dist = get_bin_dist(insert_size_mean, insert_size_stddev);

        simulator = ReadSimulator<>(output_path, ref_genome, read_size, insert_size_dist,
                                    SAM_mode, false, template_name_prefix, c_seed);
    }

    simulator.enable_SAM_writing(write_SAM);

    const auto chr_ids = get_relevant_chr_set(forest, chromosome_ids);

    const bool collect_missed = missed_SID_statistics || wide_format;

    auto result =
        simulate_seq(simulator, sequencer, forest, chr_ids, coverage,
                     with_normal_sample, purity, pre_neoplastic_in_normal,
                     true, filename_prefix, collect_missed,
                     germline_statistics, Rcpp::Rcout, c_seed, quiet);

    if (remove_output_path) {
        std::filesystem::remove_all(output_path);
    }

    using namespace Rcpp;

    auto parameters = List::create(
        _["sequencer"] = get_sequencer_data(sequencer),
        _["reference_genome"] = reference_genome, _["chromosomes"] = chromosome_ids,
        _["coverage"] = coverage, _["read_size"] = read_size,
        _["insert_size_mean"] = insert_size_mean,
        _["insert_size_stddev"] = insert_size_stddev, _["output_dir"] = output_dir,
        _["write_SAM"] = write_SAM, _["update_SAM"] = update_SAM_dir,
        _["purity"] = purity,
        _["with_normal_sample"] = with_normal_sample,
        _["filename_prefix"] = filename_prefix,
        _["template_name_prefix"] = template_name_prefix,
        _["missed_SID_statistics"] = missed_SID_statistics,
        _["germline_statistics"] = germline_statistics,
        _["wide_format"] = wide_format,
        _["seed"] = c_seed, _["quiet"] = quiet,
        _["driver_mutations"] = forest.get_driver_mutations());

    return List::create(_["mutations"] = 
                            (wide_format?get_wide_dataframe(result, missed_SID_statistics):
                                         get_long_dataframe(result)),
                        _["parameters"] = parameters);
}

Rcpp::List simulate_normal_seq(const PhylogeneticForest &forest, SEXP &sequencer,
                               SEXP &reference_genome, SEXP &chromosome_ids,
                               const double &coverage, const int &read_size,
                               const int &insert_size_mean, const int &insert_size_stddev,
                               const std::string &output_dir, const bool &write_SAM,
                               const bool &update_SAM_dir,
                               const std::string &filename_prefix,
                               const std::string &template_name_prefix,
                               const bool &missed_SID_statistics,
                               const bool &germline_statistics,
                               const bool &wide_format, const SEXP &seed,
                               const bool quiet)
{
    using namespace CLONES::Mutations::SequencingSimulations;

    const auto ref_genome = get_reference_genome(forest, reference_genome);

    std::filesystem::path output_path = output_dir;

    bool remove_output_path = false;

    if (!write_SAM) {
        remove_output_path = true;
        output_path = get_tmp_dir_path(output_dir);
    }

    ReadSimulator<>::Mode SAM_mode = ReadSimulator<>::Mode::CREATE;

    if (update_SAM_dir) {
        SAM_mode = ReadSimulator<>::Mode::UPDATE;
    }

    ReadSimulator<> simulator;
    auto c_seed = get_random_seed<int>(seed);
    if (insert_size_mean == 0) {
        simulator = ReadSimulator<>(output_path, ref_genome, read_size, SAM_mode, false,
                                    template_name_prefix, c_seed);
    } else {
        auto insert_size_dist = get_bin_dist(insert_size_mean, insert_size_stddev);

        simulator = ReadSimulator<>(output_path, ref_genome, read_size, insert_size_dist,
                                    SAM_mode, false, template_name_prefix, c_seed);
    }

    simulator.enable_SAM_writing(write_SAM);

    const auto empty_forest = forest.get_subforest_for({});

    const auto chr_ids = get_relevant_chr_set(empty_forest, chromosome_ids);

    const bool collect_missed = missed_SID_statistics || wide_format;
    auto result =
        simulate_seq(simulator, sequencer, empty_forest, chr_ids, coverage, true, 1,
                     false, true, filename_prefix, collect_missed,
                     germline_statistics, Rcpp::Rcout, c_seed, quiet);

    if (remove_output_path) {
        std::filesystem::remove_all(output_path);
    }

    using namespace Rcpp;

    auto parameters = List::create(
        _["sequencer"] = get_sequencer_data(sequencer),
        _["reference_genome"] = reference_genome, _["chromosomes"] = chromosome_ids,
        _["coverage"] = coverage, _["read_size"] = read_size,
        _["insert_size_mean"] = insert_size_mean,
        _["insert_size_stddev"] = insert_size_stddev, _["output_dir"] = output_dir,
        _["write_SAM"] = write_SAM, _["update_SAM"] = update_SAM_dir,
        _["filename_prefix"] = filename_prefix,
        _["template_name_prefix"] = template_name_prefix,
        _["missed_SID_statistics"] = missed_SID_statistics,
        _["germline_statistics"] = germline_statistics,
        _["wide_format"] = wide_format,
        _["seed"] = c_seed, _["quiet"] = quiet);

    return List::create(_["mutations"] =
                            (wide_format?get_wide_dataframe(result, missed_SID_statistics):
                                         get_long_dataframe(result)),
                        _["parameters"] = parameters);
}