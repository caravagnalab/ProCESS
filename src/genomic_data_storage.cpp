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

#include <fstream>
#include <filesystem>
#include <regex>

#include <Rcpp.h>


#include <germline.hpp>

#include <csv_reader.hpp>
#include <context_index.hpp>
#include <progress_bar.hpp>
#include <utils.hpp>


#include "genomic_data_storage.hpp"
#include "utility.hpp"


GermlineSubject::GermlineSubject()
{}

GermlineSubject::GermlineSubject(const std::string& name, const std::string& population,
                                 const std::string& super_population, const std::string& gender):
  name(name), population(population), super_population(super_population), gender(gender)
{}

Rcpp::List GermlineSubject::get_dataframe() const
{
  using namespace Rcpp;

  CharacterVector sample(1), pop(1), super_pop(1), gender(1);

  sample[0] = name;
  pop[0] = population;
  super_pop[0] = super_population;
  gender[0] = this->gender;

  return DataFrame::create(_["sample"] = sample, _["pop"] = pop,
                           _["super_pop"] = super_pop,
                           _["gender"] = gender);
}

Account::Account(std::string username, std::string password):
  username(username), password(password)
{}

GermlineStorage::GermlineStorage()
{}

GermlineStorage::GermlineStorage(const std::filesystem::path& directory):
  directory(directory)
{
  if (!std::filesystem::exists(directory)) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not exist.");
  }
  if (!std::filesystem::is_directory(directory)) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" is not a "
                             + "directory.");
  }
  if (!std::filesystem::exists(get_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"germline.csv\".");
  }
  if (!std::filesystem::exists(get_population_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"population.csv\".");
  }
  if (!std::filesystem::exists(get_population_descriptions_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"population_descr"
                             + "iptions.csv\".");
  }
  if (!std::filesystem::exists(get_alleles_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"alleles_per_chr."
                             + "csv\".");
  }
}

std::vector<GermlineSubject> GermlineStorage::get_population() const
{
  std::list<GermlineSubject> s_list;

  RACES::IO::CSVReader csv_reader(get_population_file(), true, '\t');

  for (const auto& row : csv_reader) {
    s_list.emplace_back(row.get_field(0), row.get_field(1),
                        row.get_field(2), row.get_field(3));
  }

  return {s_list.begin(), s_list.end()};
}

Rcpp::List GermlineStorage::get_population_df() const
{
  using namespace Rcpp;

  Function read_csv("read.csv");

  return read_csv(_["file"]=to_string(get_population_file()),
                  _["quote"]="", _["header"]=true, _["sep"] = "\t");
}

Rcpp::List GermlineStorage::get_population_descriptions_df() const
{
  using namespace Rcpp;

  Function read_csv("read.csv");

  return read_csv(_["file"]=to_string(get_population_descriptions_file()),
                  _["quote"]="", _["header"]=true, _["sep"] = "\t");
}

std::map<RACES::Mutations::ChromosomeId, size_t>
GermlineStorage::get_alleles_per_chromosome(const std::string& gender) const
{
  using namespace RACES::Mutations;

  std::map<ChromosomeId, size_t> alleles_per_chromosome;

  RACES::IO::CSVReader csv_reader(get_alleles_file(), true, '\t');

  const auto& header = csv_reader.get_header();
  auto found = find(header.begin(), header.end(), gender);

  if (found == header.end()) {
    throw std::runtime_error("Unknown gender " + gender + ".");
  }

  size_t index = found-header.begin();
  for (const auto& row : csv_reader) {
    auto chr_id = GenomicPosition::stochr(row.get_field(0));
    alleles_per_chromosome.insert({chr_id,
                                   std::stoul(row.get_field(index))});
  }

  return alleles_per_chromosome;
}

GermlineSubject GermlineStorage::get_subject(const std::string& subject_name) const
{
  RACES::IO::CSVReader csv_reader(get_population_file(), true, '\t');

  for (const auto& row : csv_reader) {
    if (row.get_field(0) == subject_name) {
      return {row.get_field(0), row.get_field(1),
              row.get_field(2), row.get_field(3)};
    }
  }

  throw std::runtime_error("Germline subject \"" + subject_name
                           + "\" not available.");
}


Rcpp::List GermlineStorage::get_subject_df(const std::string& subject_name) const
{
  return get_subject(subject_name).get_dataframe();
}

RACES::Mutations::GenomeMutations
GermlineStorage::build_germline(const std::string& subject_name, const bool quiet) const
{
  using namespace RACES::Mutations;

  auto bin_path = get_binary_file(subject_name);

  auto subject = get_subject(subject_name);

  auto num_of_alleles = get_alleles_per_chromosome(subject.gender);

  auto germline = GermlineMutations::load(get_file(), num_of_alleles,
                                          subject.name, Rcpp::Rcout, quiet);

  RACES::Archive::Binary::Out oarchive(bin_path);

  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  oarchive.save(germline, progress_bar, "germline");

  return germline;
}

RACES::Mutations::GenomeMutations
GermlineStorage::get_germline(const std::string& subject_name, const bool quiet) const
{
  using namespace RACES::Mutations;

  auto bin_path = get_binary_file(subject_name);

  if (!std::filesystem::exists(bin_path)) {
    return build_germline(subject_name, quiet);
  }

  RACES::Archive::Binary::In iarchive(bin_path);

  GenomeMutations germline;

  try {
    RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    iarchive.load(germline, progress_bar, "germline");
  } catch (RACES::Archive::WrongFileFormatDescr& ex) {
    raise_error(ex, "germline");
  } catch (RACES::Archive::WrongFileFormatVersion& ex) {
    raise_error(ex, "germline");
  }
  return germline;
}

bool signatures_from_COSMIC(const std::list<std::pair<std::string, std::filesystem::path>>& download_list)
{
    std::regex COSMIC_site_regex("^https://([a-zA-Z0-9_-]*).sanger.ac.uk");

    for (const auto& [src, dst]: download_list) {
      if (std::regex_search(src, COSMIC_site_regex)) {
        return true;
      }
    }

    return false;
}

inline Rcpp::Function get_function(const std::string& fun_str)
{
    using namespace Rcpp;

    // Get base R functions
    Function parse("parse");
    Function eval("eval");

    // Evaluate it in the global environment
    return eval(parse(_["text"] = fun_str));
}

Rcpp::List
Rcpp_wrap_download_list(const std::list<std::pair<std::string, std::filesystem::path>>& download_list)
{
    using namespace Rcpp;

    List wrapped_list(download_list.size());
    size_t i=0;
    for (const auto& [url, dest_filename]: download_list) {
        wrapped_list[i] = List::create(_["url"] = url,
                                       _["dest_filename"] = to_string(dest_filename));
        ++i;
    }

    return wrapped_list;
}

void download_COSMIC(const std::shared_ptr<Account>& COSMIC_account,
                     const std::list<std::pair<std::string, std::filesystem::path>>& download_list)
{
    using namespace Rcpp;

    if (COSMIC_account.get() == nullptr) {
        stop("Since April 2nd, 2025, COSMIC site (https://cancer.sanger.ac.uk/cosmic/)\n"
             "requires an account. Create an account, download SBS and ID mutation signature\n"
             "files, and pass them as parameters to `MutationEngine()` call. In alternative,\n"
             "provide COSMIC account details to `MutationEngine()` and let ProCESS download\n"
             "the signature files.\n");
    }

    const std::string r_code = "function(username, password, d_list) {\n"
                               "  if (!requireNamespace(\"rvest\", quietly = TRUE)) {\n"
                               "     stop(\"The package \\\"rvest\\\" is mandatory to download "
                               " signatures from COSMIC.\")\n  }\n"
                               "  cosmic_page <- \"https://cancer.sanger.ac.uk/cosmic/login\"\n"
                               "  cosmic_session <- rvest::session(cosmic_page)\n"
                               "  login_form <- rvest::html_form(cosmic_session)[[2]]\n"
                               "  filled_form <- rvest::html_form_set(login_form, email = username,\n"
                               "                                      pass = password)\n"
                               "  filled_form$action <- cosmic_page\n"
                               "  post_login_page <- rvest::session_submit(cosmic_session,"
                               " filled_form)\n  if (grepl(\"error while logging\",\n"
                               "            httr::content(post_login_page$response,"
                               " as=\"text\"))) {\n   stop(\"Wrong COSMIC username/password\")\n  }\n"
                               "  for (i in seq_along(d_list)) {\n"
                               "    url <- d_list[[i]]$url\n"
                               "    signature <- rvest::session_jump_to(cosmic_session, url)\n"
                               "    if (signature$response$status_code != 200) {\n"
                               "      stop(paste0(\"Cannot download file at \\\"\",\n"
                               "                  url,\"\\\".\"))\n"
                               "    }\ndest_filename <- d_list[[i]]$dest_filename\n"
                               "    writeBin(signature$response$content, dest_filename)\n}\n}\n";

    List d_list = Rcpp_wrap_download_list(download_list);

    Function rfunc = get_function(r_code);

    rfunc(COSMIC_account->get_username(), COSMIC_account->get_password(), d_list);
}

void download_file(const std::string& url, const std::filesystem::path dest_filename)
{
  if (!std::filesystem::exists(dest_filename.parent_path())) {
    std::filesystem::create_directory(dest_filename.parent_path());
  }

  using namespace Rcpp;

  // get default timeout option
  Function getOption_f("getOption");
  auto timeout = as<int>(getOption_f("timeout"));

  // raise the timeout to 1000 at least
  Function options_f("options");
  options_f(_["timeout"] = std::max(1000, timeout));

  // download the file
  Function download_f("download.file");
  download_f(_["url"] = url, _["destfile"] = to_string(dest_filename),
             _["mode"]="wb");

  // revert to the default timeout
  options_f(_["timeout"] = timeout);
}

std::filesystem::path GenomicDataStorage::download_file(const std::string& url) const
{
  const auto dest_filename = to_string(get_destination_path(url));

  ::download_file(url, dest_filename);

  return dest_filename;
}

void GenomicDataStorage::retrieve_signatures(const std::shared_ptr<Account>& COSMIC_account) const
{
    std::list<std::pair<std::string, std::filesystem::path>> download_list;

    collect_signatures_download_list<RACES::Mutations::SBSType>(download_list);
    collect_signatures_download_list<RACES::Mutations::IDType>(download_list);

    Rcpp::Rcout << "Downloading signature files..." << std::endl << std::flush;

    if (signatures_from_COSMIC(download_list)) {
        download_COSMIC(COSMIC_account, download_list);
    } else {
        for (const auto& [src, dst]: download_list) {
            ::download_file(src, dst);
        }
    }

    Rcpp::Rcout << "Signature file downloaded" << std::endl;
}

GenomicDataStorage::GenomicDataStorage(const std::shared_ptr<Account> COSMIC_account,
                                       const std::string& directory,
                                       const std::string& reference_source,
                                       const std::string& SBS_signatures_source,
                                       const std::string& indel_signatures_source,
                                       const std::string& driver_mutations_source,
                                       const std::string& passenger_CNAs_source,
                                       const std::string& germline_source):
  directory{std::filesystem::absolute(directory)}, reference_src{reference_source},
  SBS_signatures_src{SBS_signatures_source}, indel_signatures_src{indel_signatures_source},
  drivers_src{driver_mutations_source}, passenger_CNAs_src{passenger_CNAs_source},
  germline_src{germline_source}
{
  std::filesystem::create_directory(directory);

  retrieve_reference();

  retrieve_signatures(COSMIC_account);

  retrieve_drivers();
  retrieve_passenger_CNAs();
  retrieve_germline();

  const auto germline_path = get_germline_path();

  germline_storage = GermlineStorage(germline_path);
}

GenomicDataStorage::GenomicDataStorage(const std::string& directory,
                                       const std::string& reference_source,
                                       const std::string& SBS_signatures_source,
                                       const std::string& indel_signatures_source,
                                       const std::string& driver_mutations_source,
                                       const std::string& passenger_CNAs_source,
                                       const std::string& germline_source):
    GenomicDataStorage(nullptr, directory, reference_source, SBS_signatures_source,
                       indel_signatures_source, driver_mutations_source,
                       passenger_CNAs_source, germline_source)
{}

std::filesystem::path GenomicDataStorage::get_destination_path(const std::string& url) const
{
  try {
    auto name_occurrence = url.find_last_of('/')+1;

    auto filename = url.substr(name_occurrence);

    auto question_occurrence = filename.find('?');

    if (question_occurrence != std::string::npos) {
      filename = filename.substr(0, question_occurrence);
    }

    return get_directory()/filename;

  } catch (std::exception& e) {
    throw std::domain_error("\""+url+"\" is not a valid URL.");
  }
}

std::map<std::string, std::string> decompressors{
  {"gz", "gunzip"},
  {"bz2", "bunzip2"}
};

std::filesystem::path GenomicDataStorage::retrieve_reference() const
{
  if (!is_an_URL(reference_src) && !std::filesystem::exists(reference_src)) {
    throw std::runtime_error("Designed reference genome file \""
                        + reference_src + "\" does not exists.");
  }

  const auto reference_filename = get_reference_path();

  if (std::filesystem::exists(reference_filename)) {
    return reference_filename;
  }

  using namespace Rcpp;

  Rcout << "Downloading reference genome..." << std::endl << std::flush;

  auto downloaded_file = to_string(download_file(reference_src));

  Rcout << "Reference genome downloaded" << std::endl;

  auto suffix = downloaded_file.substr(downloaded_file.find_last_of('.')+1);

  if (suffix == "fa" && suffix == "fasta") {
    std::filesystem::rename(downloaded_file, reference_filename);
  } else {
    Rcout << "Decompressing reference file...";
    auto decomp_found = decompressors.find(suffix);

    if (decomp_found == decompressors.end()) {
      throw std::runtime_error("Unknown suffix \""+suffix+"\"");
    }

    Environment pkg = Environment::namespace_env("R.utils");
    Function decompress_f = pkg[decomp_found->second];

    decompress_f(_["filename"] = downloaded_file,
                 _["destname"] = to_string(reference_filename));

    Rcout << "done" << std::endl;
  }

  return reference_filename;
}

void GenomicDataStorage::retrieve_drivers() const
{
  if (std::filesystem::exists(drivers_src)) {
    return;
  }

  if (!is_an_URL(drivers_src) && !std::filesystem::exists(drivers_src)) {
    throw std::runtime_error("Designed driver mutations file \""
                             + drivers_src + "\" does not exists.");
  }

  const auto mutations_filename = get_driver_mutations_path();

  if (std::filesystem::exists(mutations_filename)) {
    return;
  }

  using namespace Rcpp;

  Rcout << "Downloading driver mutation file..." << std::endl << std::flush;

  ::download_file(drivers_src, mutations_filename);

  Rcout << "Driver mutation file downloaded" << std::endl;
}

void GenomicDataStorage::retrieve_passenger_CNAs() const
{
  if (std::filesystem::exists(passenger_CNAs_src)) {
    return;
  }

  if (!is_an_URL(passenger_CNAs_src)) {
    throw std::runtime_error("Designed passenger CNAs file \""
                             + passenger_CNAs_src
                             + "\" does not exists.");
  }

  using namespace Rcpp;

  const auto passenger_CNAs_filename = passenger_CNAs_storage_path();

  if (std::filesystem::exists(passenger_CNAs_filename)) {
    return;
  }

  Rcout << "Downloading passenger CNAs file..." << std::endl << std::flush;

  ::download_file(passenger_CNAs_src, passenger_CNAs_filename);

  Rcout << "Passenger CNAs file downloaded" << std::endl;
}

void GenomicDataStorage::retrieve_germline() const
{
  if (std::filesystem::exists(germline_src)) {
    return;
  }

  if (!is_an_URL(germline_src)) {
    throw std::runtime_error("Designed germline directory \""
                             + germline_src
                             + "\" does not exists.");
  }

  const auto germline_path = germline_storage_path();

  if (std::filesystem::exists(germline_path/"germlines.csv")) {
    return;
  }

  using namespace Rcpp;

  Rcout << "Downloading germline mutations..." << std::endl << std::flush;

  auto downloaded_file = download_file(germline_src);

  Rcout << "Germline mutations downloaded" << std::endl;

  Function untar("untar");
  untar(_["tarfile"] = to_string(downloaded_file),
        _["exdir"] = to_string(directory));
}

std::filesystem::path GenomicDataStorage::get_reference_path() const
{
    if (std::filesystem::exists(reference_src)) {
        return reference_src;
    }

    return reference_storage_path();
}

std::filesystem::path GenomicDataStorage::get_passenger_CNAs_path() const
{
    if (std::filesystem::exists(passenger_CNAs_src)) {
        return passenger_CNAs_src;
    }

    return passenger_CNAs_storage_path();
}

std::filesystem::path GenomicDataStorage::get_driver_mutations_path() const
{
    if (std::filesystem::exists(drivers_src)) {
        return drivers_src;
    }

    return driver_mutations_storage_path();
}

std::filesystem::path GenomicDataStorage::get_germline_path() const
{
    if (std::filesystem::exists(germline_src)) {
        return germline_src;
    }

    return germline_storage_path();
}

void GenomicDataStorage::save_sources() const
{
  std::ofstream of(get_directory()/"sources.csv");

  of << "reference\t" << reference_src << std::endl
     << "indel\t" << indel_signatures_src << std::endl
     << "SBS\t" << SBS_signatures_src << std::endl
     << "drivers\t" << drivers_src << std::endl
     << "passenger_CNAs\t" << passenger_CNAs_src << std::endl
     << "germline\t" << germline_src << std::endl;
}