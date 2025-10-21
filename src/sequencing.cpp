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

#include <string>

#include <Rcpp.h>

#include "seq_simulation.hpp"
#include "sequencers.hpp"

using namespace Rcpp;

RCPP_MODULE(Sequencing)
{

//' @name ErrorlessIlluminaSequencer
//' @title An error-free Illumina sequencer class
//' @description This class implements a perfect Illumina sequencers that
//'   does not commit errors.
//' @seealso `simulate_seq()`, `simulate_normal_seq()`, and
//'   `vignette("sequencing")` for usage examples
    class_<ErrorlessIlluminaSequencer>("ErrorlessIlluminaSequencer")

        .method("show", &ErrorlessIlluminaSequencer::show,
                "Show a description for the sequencer");

//' @name ErrorlessIlluminaSequencer
//' @description This method builds an error-free Illumina
//'   sequencer model.
//' @return A new error-free Illumina sequencer.
//' @examples
//' # build a sequencer model
//' sequencer <- ErrorlessIlluminaSequencer()
//' sequencer
    function("ErrorlessIlluminaSequencer", &ErrorlessIlluminaSequencer::build_sequencer,
             "Build a new error-free Illumina sequencer");

//' @name BasicIlluminaSequencer
//' @title A basic Illumina sequencer class
//' @description This class implements a basic model for Illumina sequencers.
//' @details It specifies a simulated sequencing error rate and the simulated
//'   sequencing errors will occurs according to that rate.
//' @seealso `simulate_seq()`, `simulate_normal_seq()`, and
//'   `vignette("sequencing")` for usage examples
    class_<BasicIlluminaSequencer>("BasicIlluminaSequencer")
        .method("show", &BasicIlluminaSequencer::show,
                "Show a description for the sequencer")

//' @name BasicIlluminaSequencer$error_rate
//' @title Getting the error rate
//' @description This method returns the sequencing error rate of the
//'   simulated illumina sequencer.
//' @return The sequencing error rate of the simulated sequencer.
//' @examples
//' # build a basic Illumina sequencer model whose errors occur
//' # at rate 4e-3
//' sequencer <- BasicIlluminaSequencer(4e-3)
//'
//' sequencer$error_rate
//'
//' sequencer$error_rate <- 5e-2
//'
//' sequencer$error_rate
        .property("error_rate",
                  (const double &(BasicIlluminaSequencer::*)()
                       const)(&BasicIlluminaSequencer::get_error_rate),
                  (void (BasicIlluminaSequencer::*)(const double &))(
                      &BasicIlluminaSequencer::set_error_rate),
                  "The sequencer error rate")

//' @name BasicIlluminaSequencer$random_quality_scores
//' @title Check the non-constant quality score model
//' @description This method returns `TRUE` if and only if the sequencers
//'    implements a non-constant quality score model.
//' @return `TRUE` if and only if the sequencers sequencers implements
//'    a non-constant quality score model.
//' @examples
//' # build a basic Illumina sequencer model whose quality scores are
//' # non-constant
//' sequencer <- BasicIlluminaSequencer(4e-3)
//'
//' sequencer$random_quality_scores
//'
//' sequencer$random_quality_scores <- FALSE
//'
//' sequencer$random_quality_scores
        .property("random_quality_scores",
                  (const bool &(BasicIlluminaSequencer::*)()
                       const)(&BasicIlluminaSequencer::producing_random_scores),
                  (void (BasicIlluminaSequencer::*)(const bool &))(
                      &BasicIlluminaSequencer::set_random_scores),
                  "A Boolean flag enabling non-constant quality score model");

//' @name BasicIlluminaSequencer
//' @description This method builds a basic Illumina sequencer model.
//' @param error_rate The error rate of the sequencer model.
//' @param random_quality_scores A Boolean flag to enable a basic
//'   non-constant quality score model. When it is set to `FALSE`, all
//'   the bases with no sequencing errors have the same quality score.
//'   The random quality score model increases the computation time of
//'   about 70%. (default: `TRUE`)
//' @return A basic Illumina sequencer model.
//' @examples
//' # build a sequencer model having error rate 4e-3
//' sequencer <- BasicIlluminaSequencer(error_rate=4e-3)
//' sequencer
//'
//' # build a sequencer model having error rate 4e-3 and set the seed to 5
//' sequencer <- BasicIlluminaSequencer(error_rate=4e-3, seed=5)
//' sequencer
    function("BasicIlluminaSequencer", &BasicIlluminaSequencer::build_sequencer,
             List::create(_["error_rate"], _["random_quality_scores"] = true),
             "Create a basic Illumina sequencer model");

//' @name simulate_seq
//' @title Simulating the sequencing of sampled cells
//' @description This method simulates the sequencing of the samples in a phylogenetic
//'   forest.
//' @param phylo_forest A phylogenetic forest.
//' @param sequencer The sequencer that performs the sequencing simulation
//'   (default: an `ErrorlessIlluminaSequencer`).
//' @param reference_genome The reference genome (default: NULL to use the
//'    mutation engine reference genome).
//' @param chromosomes The chromosomes that must be considered (default:
//'   `NULL`, i.e., all the reference chromosomes).
//' @param coverage The sequencing coverage (default: `10`).
//' @param read_size The read size (default: `150`).
//' @param insert_size_mean The insert size mean. Use 0 for single read
//'   sequencing and any value greater than 0 for pair read sequencing
//'   (default: `0`).
//' @param insert_size_stddev The insert size standard deviation.
//'   (default: `10`).
//' @param output_dir The SAM output directory (default:
//'   `"ProCESS_SAM"`).
//' @param write_SAM A Boolean flag to enable/disable SAM generation
//'   (default: `FALSE`).
//' @param update_SAM Update the output directory (default: `FALSE`).
//' @param purity The ratio between the number of sample tumour cell
//'   and that of all the cells, i.e., tumour and normal
//'   ones. This value must belong to the interval [0,1]
//'   (default: `1`).
//' @param with_normal_sample A Boolean flag to enable/disable the
//'   analysis of a normal sample (default: `TRUE`).
//' @param pre_neoplastic_in_normal A Boolean flag to add/remove
//'   pre-neoplastic mutations in both normal sample and normal
//'   contaminant cells (default: `FALSE`).
//' @param filename_prefix The prefix of the output SAM file name
//'   (default: `"chr_"`).
//' @param template_name_prefix The template name prefix (default:
//'   `"r"`).
//' @param missed_SID_statistics A Boolean flag to collect
//'   statistics also about the mutations that are not covered by
//'   any of the simulated reads (default: `FALSE`).
//' @param germline_statistics A Boolean flag to collect
//'   statistics also about the germinal mutations that are not
//'   covered by any of the simulated reads (default: `FALSE`).
//' @param wide_format A Boolean flag to request wide/long format
//'   for the mutation output (default: `TRUE`).
//' @param seed The random seed for the internal random generator
//'   (optional).
//' @param quiet A Boolean flag to enable/disable the progress bar
//'   (default: FALSE).
//' @return A named list of two elements: the sequencing output data
//'   frame (name "`mutations`") and the calling parameters (name
//'   "`parameters`").
//'
//'   If `wide_format` is set to `true`, the sequencing output data
//'   frame reports, for each of the observed SNVs and indels, the
//'   chromosome and the position in which it occurs (columns `chr`
//'   and `from`), the reference and the alternative sequence,
//'   the causes, and the classes of the mutation (columns `ref`,
//'   `alt`, `causes`, and `classes`, respectively). Moreover, the
//'   returned data frame contains three columns per sample: the
//'   number of reads in which the corresponding SNV occurs (column
//'   `<sample name>.NV`), the coverage of the SNV locus (column
//'   `<sample name>.DP`), and the corresponding VAF (column
//'   `<sample name>.VAF`).
//'
//'   Instead, when `wide_format` is set to `false`, the output data
//'   frame contains a row for each mutation in each sample and
//'   consists of 10 columns: `sample`, `chr`, `from`, `ref`,
//'   `alt`, `causes`, `classes`, `NV`, `DP`, and `VAF`. The column
//'   `sample` contains the name of the sample in which the
//'   mutation has been identified. The columns `chr`, from`, `ref`,
//'   `alt`, `causes`, and `classes` correspond to those of the
//'   wide_format` output. The columns `NV`, `DP`, and `VAF` 
//'   maintain the number of occurrences, the coverage, and the VAF
//'   of the mutation in cited sample.
//' @seealso `BasicIlluminaSequencer` and
//'   `ErrorlessIlluminaSequencer` as sequencer types, and
//'   `vignette("sequencing")` for usage examples
    function("simulate_seq", &simulate_seq,
             List::create(
                 _["phylo_forest"], _["sequencer"] = R_NilValue,
                 _["reference_genome"] = R_NilValue, _["chromosomes"] = R_NilValue,
                 _["coverage"] = 10, _["read_size"] = 150, _["insert_size_mean"] = 0,
                 _["insert_size_stddev"] = 10, _["output_dir"] = "ProCESS_SAM",
                 _["write_SAM"] = false, _["update_SAM"] = false, _["purity"] = 1,
                 _["with_normal_sample"] = true, _["pre_neoplastic_in_normal"] = false,
                 _["filename_prefix"] = "chr_", _["template_name_prefix"] = "r",
                 _["missed_SID_statistics"] = false, _["germline_statistics"] = false,
                 _["wide_format"] = true, _["seed"] = R_NilValue,
                 _["quiet"] = false),
             "Simulating the sequencing of the samples in a phylogenetic forest");

//' @name simulate_normal_seq
//' @title Simulating wild-type sequencing
//' @description This method simulates a wild-type sample sequencing in a
//'   phylogenetic forest. Add the cells in the wild-type sample contains
//'   the germline mutations. The forest pre-neoplastic mutations are also
//'   added to the sample by default. However, they can be avoided by
//'   using the parameter `with_pre_neoplastic`.
//' @param phylo_forest A phylogenetic forest.
//' @param sequencer The sequencer that performs the sequencing simulation
//'   (default: an `ErrorlessIlluminaSequencer`).
//' @param reference_genome The reference genome (default: NULL to use the
//'    mutation engine reference genome).
//' @param chromosomes The chromosomes that must be considered (default:
//'   `NULL`, i.e., all the reference chromosomes).
//' @param coverage The sequencing coverage (default: `10`).
//' @param read_size The read size (default: `150`).
//' @param insert_size_mean The insert size mean. Use 0 for single read
//'   sequencing and any value greater than 0 for pair read sequencing
//'   (default: `0`).
//' @param insert_size_stddev The insert size standard deviation.
//'   (default: `10`).
//' @param output_dir The SAM output directory (default:
//'   `"ProCESS_normal_SAM"`).
//' @param write_SAM A Boolean flag to enable/disable SAM generation
//'   (default: `TRUE`).
//' @param update_SAM Update the output directory (default: `FALSE`).
//' @param filename_prefix The prefix of the output SAM file name
//'   (default: `"chr_"`).
//' @param template_name_prefix The template name prefix (default:
//'   `"r"`).
//' @param missed_SID_statistics A Boolean flag to collect
//'   statistics also about the mutations that are not covered by
//'   any of the simulated reads (default: `FALSE`).
//' @param germline_statistics A Boolean flag to collect
//'   statistics also about the germinal mutations that are not
//'   covered by any of the simulated reads (default: `FALSE`).
//' @param wide_format A Boolean flag to request wide/long format
//'   for the mutation output (default: `TRUE`).
//' @param seed The random seed for the internal random generator
//'   (optional).
//' @param quiet A Boolean flag to enable/disable the progress bar
//'   (default: FALSE).
//' @return A named list of two elements: the sequencing output data
//'   frame (name "`mutations`") and the calling parameters
//'   (name "`parameters`").
//'
//'
//'   If `wide_format` is set to `true`, the sequencing output data
//'   frame reports, for each of the observed SNVs and indels, the
//'   chromosome and the position in which it occurs (columns `chr`
//'   and `from`), the reference and the alternative sequence,
//'   the causes, and the classes of the mutation (columns `ref`,
//'   `alt`, `causes`, and `classes`, respectively). Moreover, the
//'   returned data frame contains three columns: the number of
//'   reads in which the corresponding SNV occurs (column
//'   `normal.sample.NV`), the coverage of the SNV
//'   locus (column `normal.sample.DP`), and the
//'   corresponding VAF (column `normal.sample.VAF`).
//'
//'   Instead, when `wide_format` is set to `false`, the output data
//'   frame contains a row for each mutation in each sample and
//'   consists of 10 columns: `sample`, `chr`, `from`, `ref`,
//'   `alt`, `causes`, `classes`, `NV`, `DP`, and `VAF`. The column
//'   `sample` contains the name of the sample in which the
//'   mutation has been identified. The columns `chr`, from`, `ref`,
//'   `alt`, `causes`, and `classes` correspond to those of the
//'   wide_format` output. The columns `NV`, `DP`, and `VAF` 
//'   maintain the number of occurrences, the coverage, and the VAF
//'   of the mutation in cited sample.
//' @seealso `BasicIlluminaSequencer` and
//'   `ErrorlessIlluminaSequencer` as sequencer types, and
//'   `vignette("sequencing")` for usage examples
    function("simulate_normal_seq", &simulate_normal_seq,
             List::create(
                 _["phylo_forest"], _["sequencer"] = R_NilValue,
                 _["reference_genome"] = R_NilValue, _["chromosomes"] = R_NilValue,
                 _["coverage"] = 10, _["read_size"] = 150, _["insert_size_mean"] = 0,
                 _["insert_size_stddev"] = 10, _["output_dir"] = "ProCESS_normal_SAM",
                 _["write_SAM"] = true, _["update_SAM"] = false,
                 _["filename_prefix"] = "chr_", _["template_name_prefix"] = "r",
                 _["missed_SID_statistics"] = false, _["germline_statistics"] = false,
                 _["wide_format"] = true, _["seed"] = R_NilValue,
                 _["quiet"] = false),
             "Simulating the sequencing of a normal sample");

}