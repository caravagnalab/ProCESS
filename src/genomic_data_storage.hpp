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

#ifndef __PROCESS_SETUP_MUTATION_ENGINE__
#define __PROCESS_SETUP_MUTATION_ENGINE__

#include <string>

#include <Rcpp.h>

#include <utils.hpp>
#include <mutation_engine.hpp>

struct GermlineSubject
{
  std::string name;
  std::string population;
  std::string super_population;
  std::string gender;

  GermlineSubject();

  GermlineSubject(const std::string& name, const std::string& population,
                  const std::string& super_population, const std::string& gender);

  Rcpp::List get_dataframe() const;

  template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<RACES::Archive::Basic::Out, ARCHIVE>, bool> = true>
  inline void save(ARCHIVE& archive) const
  {
      archive & name
              & population
              & super_population
              & gender;
  }

  template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<RACES::Archive::Basic::In, ARCHIVE>, bool> = true>
  inline static GermlineSubject load(ARCHIVE& archive)
  {
    GermlineSubject germline_subject;

    archive & germline_subject.name
            & germline_subject.population
            & germline_subject.super_population
            & germline_subject.gender;

    return germline_subject;
  }
};

class GermlineStorage
{
  std::filesystem::path directory;

  inline std::filesystem::path get_alleles_file() const
  {
    return get_path()/std::string("alleles_per_chr.csv");
  }

  inline std::filesystem::path get_population_file() const
  {
    return get_path()/std::string("population.csv");
  }

  inline std::filesystem::path get_population_descriptions_file() const
  {
    return get_path()/std::string("population_descriptions.csv");
  }

  inline std::filesystem::path get_file() const
  {
    return get_path()/std::string("germlines.csv");
  }

  inline std::filesystem::path get_binary_file(const std::string& subject_name) const
  {
    return get_path()/("germline_" + subject_name + ".dat");
  }

  RACES::Mutations::GenomeMutations build_germline(const std::string& subject_name,
                                                   const bool quiet) const;

public:

  GermlineStorage();

  GermlineStorage(const std::filesystem::path& directory);

  inline std::filesystem::path get_path() const
  {
    return directory;
  }

  std::vector<GermlineSubject> get_population() const;

  std::map<RACES::Mutations::ChromosomeId, size_t>
  get_alleles_per_chromosome(const std::string& gender) const;

  GermlineSubject get_subject(const std::string& subject_name) const;

  RACES::Mutations::GenomeMutations get_germline(const std::string& subject_name,
                                                 const bool quiet) const;

  Rcpp::List get_subject_df(const std::string& subject_name) const;

  Rcpp::List get_population_df() const;

  Rcpp::List get_population_descriptions_df() const;
};

class GenomicDataStorage
{
  std::filesystem::path directory;
  GermlineStorage germline_storage;

  bool reference_downloaded;
  std::string reference_src;

  bool SBS_signatures_downloaded;
  std::string SBS_signatures_src;

  bool indel_signatures_downloaded;
  std::string indel_signatures_src;

  bool drivers_downloaded;
  std::string drivers_src;

  bool passenger_CNAs_downloaded;
  std::string passenger_CNAs_src;

  bool germline_downloaded;
  std::string germline_src;

  std::filesystem::path get_destination_path(const std::string& url) const;

  std::filesystem::path download_file(const std::string& url) const;

  std::filesystem::path retrieve_reference();

  static bool is_an_URL(const std::string reference)
  {
    std::set<std::string> protocols{"ftp", "http"};

    for (const auto& protocol : protocols) {
      if (reference.find(protocol)==0) {
        return true;
      }
    }

    return false;
  }

  template<typename MUTATION_TYPE>
  std::filesystem::path retrieve_file(const std::string& name, bool& downloaded)
  {
    auto source = get_signatures_path<MUTATION_TYPE>();
    if (std::filesystem::exists(source)) {
      return source;
    }

    source = get_signatures_path<MUTATION_TYPE>(false);

    downloaded = is_an_URL(source);
    if (!downloaded && !std::filesystem::exists(source)) {
      throw std::runtime_error("Designed " + name + " file \""
                               + to_string(source)
                               + "\" does not exists.");
    }

    auto filename = get_signatures_path<MUTATION_TYPE>(downloaded);
    if (std::filesystem::exists(filename)) {
      return filename;
    }

    using namespace Rcpp;

    Rcout << "Downloading " + name + " file..." << std::endl << std::flush;

    auto downloaded_file = download_file(source);

    Rcout << name << " file downloaded" << std::endl;

    std::filesystem::rename(downloaded_file, filename);

    return filename;
  }

  std::filesystem::path retrieve_indel_signatures();

  std::filesystem::path retrieve_drivers();

  std::filesystem::path retrieve_passenger_CNAs();

  std::filesystem::path retrieve_germline();
public:
  GenomicDataStorage(const std::string& directory,
                     const std::string& reference_source,
                     const std::string& SBS_signatures_source,
                     const std::string& indel_signatures_source,
                     const std::string& driver_mutations_source,
                     const std::string& passenger_CNAs_source,
                     const std::string& germline_source);

  inline std::filesystem::path get_directory() const
  {
    return directory;
  }

  std::filesystem::path get_reference_path() const;

  template<typename MUTATION_TYPE,
           std::enable_if_t<std::is_base_of_v<RACES::Mutations::MutationType, MUTATION_TYPE>, bool> = true>
  std::string get_signatures_path(const bool& downloaded) const
  {
    if constexpr(std::is_base_of_v<MUTATION_TYPE, RACES::Mutations::SBSType>) {
      if (downloaded) {
        return to_string(directory/std::string("SBS_signatures.txt"));
      } else {
        return SBS_signatures_src;
      }
    }

    if constexpr(std::is_base_of_v<MUTATION_TYPE, RACES::Mutations::IDType>) {
      if (downloaded) {
        return to_string(directory/std::string("indel_signatures.txt"));
      } else {
        return indel_signatures_src;
      }
    }

    throw std::runtime_error("Unsupported mutation type.");
  }

  template<typename MUTATION_TYPE,
           std::enable_if_t<std::is_base_of_v<RACES::Mutations::MutationType, MUTATION_TYPE>, bool> = true>
  std::string get_signatures_path() const
  {
    if constexpr(std::is_base_of_v<MUTATION_TYPE, RACES::Mutations::SBSType>) {
      return get_signatures_path<MUTATION_TYPE>(SBS_signatures_downloaded);
    }

    if constexpr(std::is_base_of_v<MUTATION_TYPE, RACES::Mutations::IDType>) {
      return get_signatures_path<MUTATION_TYPE>(indel_signatures_downloaded);
    }

    throw std::runtime_error("Unsupported mutation type.");
  }

  std::filesystem::path get_driver_mutations_path() const;

  std::filesystem::path get_passenger_CNAs_path() const;

  inline const GermlineStorage& get_germline_storage() const
  {
    return germline_storage;
  }

  void save_sources() const;
};

#endif // __PROCESS_SETUP_MUTATION_ENGINE__
