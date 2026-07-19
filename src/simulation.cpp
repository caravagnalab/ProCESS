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

#include <algorithm>
#include <filesystem>

#include <ending_conditions.hpp>

#include "simulation.hpp"

#include "utility.hpp"

const std::map<CLONES::Mutants::CellEventType, std::string> event2name = {
    {CLONES::Mutants::CellEventType::DEATH, "death"},
    {CLONES::Mutants::CellEventType::DUPLICATION, "duplication"},
    {CLONES::Mutants::CellEventType::MUTATION, "mutation"},
    {CLONES::Mutants::CellEventType::DUP_AND_EPI_SWITCH, "switch"}
};

template<typename KEY, typename VALUE>
std::map<VALUE, KEY> invert_map(const std::map<KEY, VALUE>& original)
{
    std::map<VALUE, KEY> inverted;

    for (const auto& [key, value]: original) {
        inverted.emplace(value, key);
    }

    return inverted;
}

const std::map<std::string, CLONES::Mutants::CellEventType> name2event{invert_map(event2name)};

inline void
handle_unknown_event(const std::string &event,
                     const std::map<std::string, CLONES::Mutants::CellEventType>& name2event)
{
    std::ostringstream oss;

    oss << "Event descriptor \"" << event << "\" is not supported. " << std::endl
        << "Supported descriptors are ";

    size_t i{0};
    for (const auto &[name, type] : name2event) {
        if (type != CLONES::Mutants::CellEventType::MUTATION) {
            if (i > 0) {
                if (name2event.size() != 2) {
                    oss << ",";
                }
                oss << " ";
            }

            if ((++i) == name2event.size()) {
                oss << "and ";
            }

            oss << "\"" << name << "\"";
        }
    }

    oss << ".";

    Rcpp::stop(oss.str());
}

const CLONES::Mutants::CellEventType&
get_event_id(const std::string& event_name,
             const std::map<std::string, CLONES::Mutants::CellEventType>& name2event)
{
    auto event_it = name2event.find(event_name);

    if (event_it == name2event.end()) {
        handle_unknown_event(event_name, name2event);
    }

    return event_it->second;
}

inline const CLONES::Mutants::CellEventType& get_event_id(const std::string& event_name)
{
    return get_event_id(event_name, name2event);
}

const std::string&
get_event_name(const CLONES::Mutants::CellEventType& event_id,
               const std::map<CLONES::Mutants::CellEventType, std::string>& event2name)
{
    auto event_it = event2name.find(event_id);

    if (event_it == event2name.end()) {
        Rcpp::stop("Unknown event identifier "
                   + std::to_string(static_cast<uint>(event_id)) + ".");
    }

    return event_it->second;
}

inline const std::string& get_event_name(const CLONES::Mutants::CellEventType& event_id)
{
    return get_event_name(event_id, event2name);
}

template <typename SIMULATION_TEST> struct RTest : public SIMULATION_TEST
{
    size_t counter;

    template <typename... Args>
    explicit RTest(Args... args) : SIMULATION_TEST(args...), counter{0}
    {}

    bool operator()(const CLONES::Mutants::Evolutions::TissueSimulation &simulation)
    {
        if (++counter >= 10000) {
            counter = 0;
            try {
                Rcpp::checkUserInterrupt();
            } catch (...) {
                return true;
            }
        }

        using namespace CLONES::Mutants::Evolutions;

        return SIMULATION_TEST::operator()(simulation);
    }
};

size_t count_events(const CLONES::Mutants::Evolutions::SpeciesStatistics &statistics,
                    const CLONES::Mutants::CellEventType &event)
{
    switch (event) {
    case CLONES::Mutants::CellEventType::DEATH:
        return statistics.killed_cells;
    case CLONES::Mutants::CellEventType::DUPLICATION:
        return statistics.num_duplications;
    case CLONES::Mutants::CellEventType::DUP_AND_EPI_SWITCH:
        return statistics.num_of_epigenetic_events();
    default:
        Rcpp::stop("get_counts: unsupported event");
    }
}

std::set<CLONES::Mutants::SpeciesId>
get_species_ids_from_mutant_name(const CLONES::Mutants::Evolutions::Tissue &tissue,
                                 const std::set<std::string> &mutant_name)
{
    std::set<CLONES::Mutants::SpeciesId> species_ids;

    for (const auto &species : tissue) {
        if (mutant_name.count(species.get_mutant_name()) > 0) {
            species_ids.insert(species.get_id());
        }
    }

    return species_ids;
}

CLONES::Mutants::Evolutions::PositionInTissue get_position_in_tissue(
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &position)
{
    if (position.size() == 2) {
        return {position[0], position[1]};
    }

    Rcpp::stop("ProCESS supports only the 2D space.");
}

CLONES::Mutants::RectangleSet
TissueSimulation::get_rectangle(const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
                                const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner)
{
    auto l_position = get_position_in_tissue(lower_corner);
    auto u_position = get_position_in_tissue(upper_corner);

    return {l_position, u_position};
}

size_t count_driver_mutated_cells(
    const CLONES::Mutants::Evolutions::Tissue &tissue,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner,
    const std::set<CLONES::Mutants::SpeciesId> &species_filter,
    const std::set<std::string> &epigenetic_filter)
{
    using namespace CLONES::Mutants;
    using namespace CLONES::Mutants::Evolutions;

    if (lower_corner.size() != upper_corner.size()) {
        Rcpp::stop("lower_corner and upper_corner must have the same size.");
    }

    auto lower_it = lower_corner.begin();
    auto upper_it = upper_corner.begin();
    for (; lower_it != lower_corner.end(); ++lower_it, ++upper_it) {
        if (*lower_it > *upper_it) {
            return 0;
        }
    }

    size_t total{0};
    for (auto x = lower_corner[0]; x <= upper_corner[0]; ++x) {
        for (auto y = lower_corner[1]; y <= upper_corner[1]; ++y) {
            auto cell_proxy = tissue({x, y});
            if (!cell_proxy.is_wild_type()) {
                const CellInTissue &cell = cell_proxy;

                if (species_filter.count(cell.get_species_id()) > 0) {

                    const auto &species = tissue.get_species(cell.get_species_id());
                    if (tissue.get_epigenetic_state_names().size() == 0
                        || epigenetic_filter.count(species.get_epistate_name()) > 0) {
                        ++total;
                    }
                }
            }
        }
    }

    return total;
}

std::vector<CLONES::Mutants::Evolutions::Direction>
TissueSimulation::get_possible_directions()
{
    namespace RS = CLONES::Mutants::Evolutions;

    std::vector<RS::Direction> directions;
    for (const auto &x_move :
         {RS::Direction::X_UP, RS::Direction::X_DOWN, RS::Direction::X_NULL}) {
        for (const auto &y_move :
             {RS::Direction::Y_UP, RS::Direction::Y_DOWN, RS::Direction::Y_NULL}) {
            directions.push_back(x_move | y_move);
        }
    }

    // remove null move
    directions.pop_back();

    return directions;
}

PlainChooser::PlainChooser(
    const std::shared_ptr<CLONES::Mutants::Evolutions::TissueSimulation> &sim_ptr,
    const std::string &mutant_name)
    : sim_ptr(sim_ptr), mutant_name(mutant_name)
{}

RectangularChooser::RectangularChooser(
    const std::shared_ptr<CLONES::Mutants::Evolutions::TissueSimulation> &sim_ptr,
    const std::string &mutant_name,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner)
    : PlainChooser(sim_ptr, mutant_name),
      rectangle(TissueSimulation::get_rectangle(lower_corner, upper_corner))
{}

bool TissueSimulation::has_names(const Rcpp::List &list,
                                 std::vector<std::string> aimed_names)
{
    if (aimed_names.size() != static_cast<size_t>(list.size())) {
        return false;
    }

    for (const std::string &name : aimed_names) {
        if (!list.containsElementNamed(name.c_str())) {
            return false;
        }
    }

    return true;
}

bool TissueSimulation::has_names_in(const Rcpp::List &list,
                                    std::set<std::string> aimed_names)
{
    using namespace Rcpp;

    if (aimed_names.size() < static_cast<size_t>(list.size())) {
        return false;
    }

    CharacterVector names = wrap(list.names());

    for (size_t i = 0; i < static_cast<size_t>(names.size()); ++i) {
        if (aimed_names.count(as<std::string>(names[i])) == 0) {
            return false;
        }
    }

    return true;
}

Rcpp::List TissueSimulation::get_cells(const std::string &sample_name) const
{
    if (!already_collected_sample(sample_name)) {
        Rcpp::stop("Unknown sample \"" + sample_name + "\".");
    }

    const auto tissue = load_tissue(pre_sample_tissue_image_path(sample_name));

    return get_cells(tissue);
}

inline Rcpp::String encode_epigenetic_state(const CLONES::Mutants::Evolutions::Species& species)
{
    using namespace Rcpp;

    if (species.get_epistate_name() == "") {
        return NA_STRING;
    }

    return species.get_epistate_name();
}

Rcpp::List TissueSimulation::get_cells(
    const CLONES::Mutants::Evolutions::Tissue &tissue,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner,
    const std::set<CLONES::Mutants::SpeciesId> &species_filter,
    const std::set<std::string> &epigenetic_filter)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    namespace RS = CLONES::Mutants::Evolutions;

    if (lower_corner.size() != 2) {
        Rcpp::stop("The lower corner must be a vector having size 2");
    }

    if (upper_corner.size() != 2) {
        Rcpp::stop("The upper corner must be a vector having size 2");
    }

    size_t num_of_rows = count_driver_mutated_cells(tissue, lower_corner, upper_corner,
                                                    species_filter, epigenetic_filter);

    IntegerVector ids(num_of_rows);
    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
    IntegerVector x_pos(num_of_rows), y_pos(num_of_rows);

    size_t i{0};
    for (auto x = lower_corner[0]; x <= upper_corner[0]; ++x) {
        for (auto y = lower_corner[1]; y <= upper_corner[1]; ++y) {
            auto cell_proxy = tissue({x, y});
            if (!cell_proxy.is_wild_type()) {

                const RS::CellInTissue &cell = cell_proxy;

                const auto &species = tissue.get_species(cell.get_species_id());

                if (species_filter.count(cell.get_species_id()) > 0
                        && (tissue.get_epigenetic_state_names().size() == 0
                            || epigenetic_filter.count(species.get_epistate_name()) > 0)) {

                    ids[i] = cell.get_id();
                    mutant_names[i] = species.get_mutant_name();

                    epi_states[i] = encode_epigenetic_state(species);

                    x_pos[i] = x;
                    y_pos[i] = y;

                    ++i;
                }
            }
        }
    }

    if (tissue.get_epigenetic_state_names().size() == 0) {
        return DataFrame::create(_["cell_id"] = ids, _["mutant"] = mutant_names,
                                 _["position_x"] = x_pos, _["position_y"] = y_pos);
    }

    return DataFrame::create(_["cell_id"] = ids, _["mutant"] = mutant_names,
                             _["epistate"] = epi_states, _["position_x"] = x_pos,
                             _["position_y"] = y_pos);

}

Rcpp::List
TissueSimulation::wrap_a_cell(const CLONES::Mutants::Evolutions::CellInTissue &cell) const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    const auto &species = sim_ptr->tissue().get_species(cell.get_species_id());

    const auto mutant_name = species.get_mutant_name();
    const auto epistate = species.get_epistate_name();

    if (epistate != "") {
        return DataFrame::create(_["cell_id"] = cell.get_id(), _["mutant"] = mutant_name,
                                _["epistate"] = epistate, _["position_x"] = cell.x,
                                _["position_y"] = cell.y);
    }

    return DataFrame::create(_["cell_id"] = cell.get_id(), _["mutant"] = mutant_name,
                             _["position_x"] = cell.x, _["position_y"] = cell.y);
}

TissueSimulation TissueSimulation::load(const std::string &directory_name)
{
    using namespace CLONES::Mutants::Evolutions;

    TissueSimulation simulation;

    simulation.save_snapshots = true;
    simulation.name = directory_name;

    if (!std::filesystem::exists(directory_name)) {
        Rcpp::stop("The directory \"" + directory_name + "\" does not exist.");
    }

    if (!std::filesystem::is_directory(directory_name)) {
        Rcpp::stop("\"" + directory_name + "\" is not a directory.");
    }

    auto snapshot_path = BinaryLogger::find_last_snapshot_in(directory_name);

    CLONES::Archive::Binary::In archive(snapshot_path);

    try {
        archive &*(simulation.sim_ptr);
    } catch (const CLONES::Archive::WrongFileFormatDescr &ex) {
        raise_error(ex, "tissue simulation");
    } catch (const CLONES::Archive::WrongFileFormatVersion &ex) {
        raise_error(ex, "tissue simulation");
    }

    auto ruh_path =
        std::filesystem::path(directory_name) / get_rates_update_history_file_name();

    if (std::filesystem::exists(ruh_path)) {
        CLONES::Archive::Binary::In ruh_archive(ruh_path);

        ruh_archive & simulation.rate_update_history;
    } else {
        Rcpp::warning("The rates update history file is missing.");
    }
    return simulation;
}

std::string get_time_string()
{
    std::time_t time;
    std::tm *info;
    char buffer[81];

    std::time(&time);
    info = std::localtime(&time);

    std::strftime(buffer, 80, "%Y%m%d-%H%M%S", info);

    return buffer;
}

inline std::string get_default_name() { return "ProCESS_" + get_time_string(); }

void TissueSimulation::init_history_rate_updates()
{
    const auto dir_path = std::filesystem::absolute(sim_ptr->get_logger().get_directory());
    if (!std::filesystem::exists(dir_path)) {
        std::filesystem::create_directory(sim_ptr->get_logger().get_directory());
    }

    // save changes
    CLONES::Archive::Binary::Out ruh_archive(get_rates_update_history_path());

    ruh_archive & rate_update_history;
}

void TissueSimulation::init(const SEXP &sexp)
{
    using namespace Rcpp;

    namespace RS = CLONES::Mutants::Evolutions;

    switch (TYPEOF(sexp)) {
        case INTSXP:
        case REALSXP:
        {
            int seed = as<int>(sexp);
            name = get_default_name();

            if (save_snapshots) {
                sim_ptr = std::make_shared<RS::TissueSimulation>(name, seed);
            } else {
                sim_ptr = std::make_shared<RS::TissueSimulation>(get_tmp_dir_path(), seed);
            }
            break;
        }
        case STRSXP:
        {
            name = as<std::string>(sexp);

            auto seed = get_random_seed<int>(R_NilValue);
            if (save_snapshots) {
                sim_ptr = std::make_shared<RS::TissueSimulation>(name, seed);
            } else {
                sim_ptr = std::make_shared<RS::TissueSimulation>(get_tmp_dir_path(), seed);
            }
            break;
        }
        default:
        {
            std::ostringstream oss;

            oss << "Invalid type for the first parameter: " << type2name(sexp);

            Rcpp::stop(oss.str());
        }
    }

    init_history_rate_updates();
}

TissueSimulation::TissueSimulation()
    : sim_ptr{std::make_shared<CLONES::Mutants::Evolutions::TissueSimulation>(
          get_tmp_dir_path(), get_random_seed<int>(R_NilValue))},
      name{get_default_name()}, save_snapshots{false}
{
    TissueSimulation::cell_event_names_inv = invert_map(CLONES::Mutants::cell_event_names);

    init_history_rate_updates();
}

TissueSimulation::TissueSimulation(const SEXP &sexp) : save_snapshots{false}
{
    TissueSimulation::cell_event_names_inv = invert_map(CLONES::Mutants::cell_event_names);

    using namespace Rcpp;

    if (TYPEOF(sexp) == LGLSXP) {
        save_snapshots = as<bool>(sexp);
        name = get_default_name();

        auto seed = get_random_seed<int>(R_NilValue);

        std::filesystem::path sim_path;
        if (save_snapshots) {
            sim_path = name;
        } else {
            sim_path = get_tmp_dir_path();
        }
        sim_ptr = std::make_shared<CLONES::Mutants::Evolutions::TissueSimulation>(
                sim_path, seed);

        init_history_rate_updates();

        return;
    }

    init(sexp);
}

TissueSimulation::TissueSimulation(const SEXP &first_param, const SEXP &second_param)
    : save_snapshots{false}
{
    TissueSimulation::cell_event_names_inv = invert_map(CLONES::Mutants::cell_event_names);

    using namespace Rcpp;

    if (TYPEOF(second_param) == LGLSXP) {
        save_snapshots = as<bool>(second_param);

        init(first_param);

        return;
    }

    if (TYPEOF(first_param) != STRSXP) {
        std::ostringstream oss;

        oss << "Invalid type for the parameter 1: " << type2name(first_param)
            << ". If the last parameter is not a Boolean value (save on disk"
            << " parameter), it must be a string (the name of the simulation).";

        Rcpp::stop(oss.str());
    }

    if (TYPEOF(second_param) != INTSXP && TYPEOF(second_param) != REALSXP) {
        std::ostringstream oss;

        oss << "Invalid type for the parameter 2: " << type2name(second_param)
            << ". If the last parameter is not a Boolean value (save on disk"
            << " parameter), it must be an integer value (the random seed).";

        Rcpp::stop(oss.str());
    }

    name = as<std::string>(first_param);
    int seed = as<int>(second_param);

    sim_ptr = std::make_shared<CLONES::Mutants::Evolutions::TissueSimulation>(name, seed);

    init_history_rate_updates();
}

TissueSimulation::TissueSimulation(const std::string &simulation_name, const SEXP &seed,
                                   const bool &save_snapshots)
    : TissueSimulation{simulation_name, get_random_seed<int>(seed), save_snapshots}
{}

TissueSimulation::TissueSimulation(const std::string &simulation_name, const int seed,
                                   const bool &save_snapshots)
    : name{simulation_name}, save_snapshots{save_snapshots}
{
    TissueSimulation::cell_event_names_inv = invert_map(CLONES::Mutants::cell_event_names);

    std::filesystem::path sim_path;
    if (save_snapshots) {
        sim_path = simulation_name;
    } else {
        sim_path = get_tmp_dir_path();
    }

    sim_ptr = std::make_shared<CLONES::Mutants::Evolutions::TissueSimulation>(
                sim_path, seed);

    init_history_rate_updates();
}

bool has_column(const Rcpp::DataFrame &df, const std::string& col_name) {
    Rcpp::CharacterVector colnames = df.attr("names");

    const auto it = std::find(colnames.begin(), colnames.end(), col_name);

    return it != colnames.end();
}

std::set<std::string> collect_mutants(const Rcpp::DataFrame& df_rates)
{
    using namespace Rcpp;

    std::set<std::string> mutants;

    if (!has_column(df_rates, "mutant")) {
        Rcpp::stop("The data frame misses the column \"mutant\".");
    }

    CharacterVector mutant_col = df_rates["mutant"];

    for (auto i=0; i<mutant_col.size(); ++i) {
        if (!CharacterVector::is_na(mutant_col[i])) {
            mutants.insert(as<std::string>(mutant_col[i]));
        }
    }

    return mutants;
}

std::set<std::string> collect_epigenetic_states(const Rcpp::DataFrame& df_rates)
{
    using namespace Rcpp;

    std::set<std::string> epi_states;

    for (const auto col_name : {"epistate", "first.child.epistate"}) {
        if (has_column(df_rates, col_name)) {
            CharacterVector epistates = df_rates[col_name];

            for (auto i=0; i<epistates.size(); ++i) {
                if (!CharacterVector::is_na(epistates[i])) {
                    epi_states.insert(as<std::string>(epistates[i]));
                }
            }
        }
    }

    return epi_states;
}

std::list<std::string> collect_strings(const Rcpp::List& string_list)
{
    using namespace Rcpp;

    std::list<std::string> strings;
    for (auto i=0; i<string_list.size(); ++i) {
        if (!is<String>(string_list[i])) {
            throw CLONES::Error<std::domain_error>("The " + ordinal_suffix(i+1)
                                                   + " element is not a string.");
        }

        strings.push_back(as<std::string>(string_list[i]));
    }

    return strings;
}

TissueSimulation TissueSimulation::build_simulation(const SEXP &simulation_name,
                                                    const SEXP &width, const SEXP &height,
                                                    const SEXP &save_snapshots,
                                                    const SEXP &rates,
                                                    const SEXP &epistates,
                                                    const SEXP &seed)
{
    std::string c_name;
    auto c_width = FromSEXP::get<size_t>(width, "parameter \"width\"", "positive natural value");
    auto c_height = FromSEXP::get<size_t>(height, "parameter \"height\"", "positive natural value");
    auto c_save = FromSEXP::get<bool>(save_snapshots, "parameter \"save_snapshots\"", "Boolean value");
    auto c_seed = get_random_seed<int>(seed);

    if (TYPEOF(simulation_name) == NILSXP) {
        c_name = get_default_name();
    } else {
        c_name = FromSEXP::get<std::string>(simulation_name, "parameter \"name\"", "string");
    }

    TissueSimulation sim(c_name, c_seed, c_save);

    sim.update_tissue(c_width, c_height);

    std::set<std::string> C_epistates;

    if (TYPEOF(epistates) != NILSXP) {
        try {
            for (const auto& epistate : collect_strings(epistates)) {
                if (C_epistates.count(epistate)>0) {
                    Rcpp::warning("Epistate \"" + epistate
                                  +"\" is listed more than once.");
                }
                C_epistates.insert(epistate);
            }
        } catch(const std::exception& ex) {
            std::string msg = "\"epigenetic_states\" must be a list of "
                              "epigenetic states.";
            if (strlen(ex.what())>0) {
                msg += " " + std::string(ex.what());
            }
            Rcpp::stop(msg);
        }
    }

    if (TYPEOF(rates) != NILSXP) {
        using namespace Rcpp;
        if (!is<DataFrame>(rates)) {
            stop("\"rates\" must be a data frame.");
        }

        DataFrame df_rates = as<DataFrame>(rates);

        for (const auto& mutant : collect_mutants(df_rates)) {
            sim.add_mutant(mutant);
        }

        auto df_epistates = collect_epigenetic_states(df_rates);

        for (const auto& epi_state : collect_epigenetic_states(df_rates)) {
            if (C_epistates.size()>0 && C_epistates.count(epi_state)>0) {
                warning("\"" + epi_state + "\" is not among the specified "
                        "epigenetic states.");
            } else {
                C_epistates.erase(epi_state);
            }
            sim.add_epistate(epi_state);
        }

        sim.set_rates(df_rates);
    }

    for (const auto& epi_state : C_epistates) {
        sim.add_epistate(epi_state);
    }

    return sim;
}

TissueSimulation::~TissueSimulation()
{
    if (sim_ptr.use_count() == 1 && !save_snapshots) {
        auto dir = sim_ptr->get_logger().get_directory();

        sim_ptr = nullptr;

        std::filesystem::remove_all(dir);
    }
}

std::set<std::string> get_epistate_names(const Rcpp::List& list, const char* what)
{
    if (!list.hasAttribute("names")) {
        std::ostringstream oss;

        oss << "The list of the " << what << " rates should be a named list whose "
            << "names are the epigenetic state names.";

        Rcpp::stop(oss.str());
    }

    const auto name_list = Rcpp::as<std::list<std::string>>(list.attr("names"));

    std::set<std::string> epistate_names{name_list.begin(), name_list.end()};

    if (epistate_names.count("")>0) {
        std::ostringstream oss;

        oss << "The list of the " << what
            << " rates contains \"\" among its names. Epigenetic states cannot "
            << "have an empty name.";

        Rcpp::stop(oss.str());
    }

    return epistate_names;
}

std::set<std::string> get_epistate_names(const Rcpp::DataFrame& df)
{
    std::set<std::string> epistate_names;

    for (const char* what: {"source", "destination"})
    {
        if (!df.containsElementNamed(what)) {
            std::ostringstream oss;

            oss << "The epigenetic switches dataframe must contain the column \""
                << what << "\".";
            Rcpp::stop(oss.str());
        }

        const auto names = Rcpp::as<std::vector<std::string>>(df[what]);

        for (size_t i=0; i<names.size(); ++i) {
            if (names[i].size()) {
                std::ostringstream oss;

                oss << "The epigenetic state name \"\" is not allowed "
                    << "(see the " << i+1 << "-" << ordinal_suffix(i+1)
                    << " row of the epigenetic switch dataframe).";

                Rcpp::stop(oss.str());
            }
        }

        epistate_names.insert(names.begin(), names.end());
    }

    return epistate_names;
}

void TissueSimulation::add_mutant(const std::string &mutant_name)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    if (mutant_name == "Wild-type") {
        Rcpp::stop("\"Wild-type\" is a reserved mutant name.");
    }

    try {
        sim_ptr->add_mutant(mutant_name);
    } catch (const std::exception& ex) {
        Rcpp::stop(ex.what());
    }
}

void TissueSimulation::add_mutant(const SEXP &mutant_name)
{
    const auto C_mutant_name = FromSEXP::get<std::string>(mutant_name, "parameter", "string");

    add_mutant(C_mutant_name);
}

void TissueSimulation::add_mutant(const std::string &mutant_name, const Rcpp::List& rates)
{
    add_mutant(mutant_name);

    set_rates(mutant_name, rates);
}

Rcpp::List char_vector_to_list(Rcpp::NumericVector orig)
{
    Rcpp::List L(orig.size());
    for(auto i = 0; i < orig.size(); ++i) {
        L[i] = orig[i];
    }

    L.names() = orig.names();

    return L;
}

void TissueSimulation::add_mutant(const SEXP &mutant_name, const SEXP& rates)
{
    const auto C_mutant_name = FromSEXP::get<std::string>(mutant_name, "1st parameter", "string");

    bool valid{false};

    Rcpp::List C_rates;
    if (Rcpp::is<Rcpp::NumericVector>(rates)) {
        auto nv = Rcpp::as<Rcpp::NumericVector>(rates);
        if (nv.hasAttribute("names")) {
            valid = true;
            C_rates = char_vector_to_list(nv);
        }

        Rcpp::CharacterVector cv = C_rates.names();

        for (auto i = 0; i < cv.size(); ++i) {
            if (cv[i] == "") {
                Rcpp::stop("If the second parameter is a numeric vector,"
                           "all its values must have a name.");
            }
        }
    } else {
        if (Rcpp::is<Rcpp::List>(rates)) {
            valid = true;
            C_rates = Rcpp::as<Rcpp::List>(rates);
        }
    }

    if (!valid) {
        Rcpp::stop("The 2nd parameter must be either a named list "
                   "or a named vector.");
    }

    add_mutant(C_mutant_name);

    set_rates(C_mutant_name, C_rates);
}

void TissueSimulation::add_epistate(const std::string &epistate_name)
{
    using namespace Rcpp;
    using namespace CLONES;
    using namespace CLONES::Mutants;

    if (epistate_name == "Wild-type") {
        Rcpp::stop("\"Wild-type\" is a reserved epigenetic name.");
    }

    try {
        sim_ptr->add_epigenetic_state(epistate_name);
    } catch (const std::exception& ex) {
        Rcpp::stop(ex.what());
    }
}

Rcpp::List TissueSimulation::get_rates() const
{
    using namespace Rcpp;

    size_t num_of_rows{0};
    for (const auto &species : sim_ptr->tissue()) {
        for (const auto &[event_id, event_rates] : species.get_rates()) {
            num_of_rows += event_rates.size();
        }
    }

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows),
                    event_names(num_of_rows), dst_epi_states(num_of_rows);
    NumericVector rates(num_of_rows);

    using namespace CLONES::Mutants;

    const auto& tissue = sim_ptr->tissue();
    size_t i{0};
    for (const auto &species : tissue) {
        for (const auto &[event_id, event_rates] : species.get_rates()) {
            for (const auto &[dst_id, rate] : event_rates) {
                mutant_names[i] = species.get_mutant_name();
                epi_states[i] = encode_epigenetic_state(species);
                event_names[i] = get_event_name(event_id);
                const auto &dst_species = tissue.get_species(dst_id);
                if (event_id==CellEventType::DEATH) {
                    dst_epi_states[i] = NA_STRING;
                } else {
                    dst_epi_states[i] = encode_epigenetic_state(dst_species);
                }
                rates[i] = rate;

                ++i;
            }
        }
    }

    if (sim_ptr->get_epigenetic_state_names().size()==0) {
        return DataFrame::create(_["mutant"] = mutant_names,
                                 _["event"] = event_names,
                                 _["rate"] = rates);
    }

    return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states,
                             _["event"] = event_names,
                             _["first child epistate"] = dst_epi_states,
                             _["rate"] = rates);
}

Rcpp::List TissueSimulation::get_rates(const SEXP& complete) const
{
    using namespace Rcpp;

    if (!(Rf_isLogical(complete) && Rf_length(complete) == 1)) {
        stop("The parameter must be a Boolean value.");
    }

    const bool C_complete = as<bool>(complete);

    if (!C_complete) {
        return get_rates();
    }

    const auto epistates = sim_ptr->get_epigenetic_state_names();
    const size_t num_epistates{epistates.size()};

    const size_t num_mutants{sim_ptr->get_mutant_names().size()};
    const size_t num_of_rows = (num_epistates>0?num_mutants*num_epistates*(num_epistates+1):
                                    num_mutants*2);

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows),
                    event_names(num_of_rows), dst_epi_states(num_of_rows);
    NumericVector rates(num_of_rows);

    using namespace CLONES::Mutants;

    const auto& tissue = sim_ptr->tissue();
    size_t i{0};
    for (const auto &species : tissue) {
        const auto& species_rates = species.get_rates();

        const auto mutant = species.get_mutant_name();
        const auto epistate = encode_epigenetic_state(species);

        // DUPLICATION AND DEATH
        for (const auto event_id : {CellEventType::DUPLICATION,
                                    CellEventType::DEATH}) {
            mutant_names[i] = mutant;
            event_names[i] = get_event_name(event_id);
            if (num_epistates>0) {
                epi_states[i] = epistate;
                if (event_id == CellEventType::DEATH) {
                    dst_epi_states[i] = NA_STRING;
                } else {
                    dst_epi_states[i] = epistate;
                }
            }
            auto rate_it = species_rates.find(event_id);
            if (rate_it == species_rates.end()) {
                rates[i] = 0;
            } else {
                rates[i] = rate_it->second.begin()->second;
            }
            ++i;
        }

        if (num_epistates>0) {
            // SWITCH
            auto rate_it = species_rates.find(CellEventType::DUP_AND_EPI_SWITCH);
            for (const auto& mutant_species : tissue.get_mutant_view(mutant)) {
                if (species.get_id() != mutant_species.get_id()) {
                    mutant_names[i] = mutant;
                    epi_states[i] = epistate;
                    event_names[i] = get_event_name(CellEventType::DUP_AND_EPI_SWITCH);
                    dst_epi_states[i] = mutant_species.get_epistate_name();
                    if (rate_it == species_rates.end()) {
                        rates[i] = 0;
                    } else {
                        const auto& dst_rates = rate_it->second;
                        auto dst_it = dst_rates.find(mutant_species.get_id());
                        if (dst_it == dst_rates.end()) {
                            rates[i] = 0;
                        } else {
                            rates[i] = dst_it->second;
                        }
                    }
                    ++i;
                }
            }
        }
    }

    if (num_epistates==0) {
        return DataFrame::create(_["mutant"] = mutant_names,
                                 _["event"] = event_names,
                                 _["rate"] = rates);
    }

    return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states,
                             _["event"] = event_names,
                             _["first child epistate"] = dst_epi_states,
                             _["rate"] = rates);
}

template<typename RETURNED_TYPE, typename CONTAINER>
RETURNED_TYPE fill_vector(const CONTAINER& container)
{
    RETURNED_TYPE result(container.size());

    size_t i{0};
    for (const auto &value : container) {
        result[i] = value;

        ++i;
    }

    return result;
}

Rcpp::List TissueSimulation::get_mutant_names() const
{
    using namespace Rcpp;

    const auto mutant_names = sim_ptr->get_mutant_names();

    auto names = fill_vector<CharacterVector>(mutant_names);

    return DataFrame::create(_["mutant"] = names);
}

Rcpp::List TissueSimulation::get_epigenetic_state_names() const
{
    using namespace Rcpp;

    const auto epistate_names = sim_ptr->get_epigenetic_state_names();

    auto names = fill_vector<CharacterVector>(epistate_names);

    return DataFrame::create(_["epistate"] = names);
}

Rcpp::List TissueSimulation::get_rates_update_history() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    size_t num_of_rows{0};
    for (const auto &[time, species_rate_updates] : rate_update_history) {
        for (const auto &[species_id, event_rate_updates] : species_rate_updates) {
            for (const auto &[event_name, dst_rate] : event_rate_updates) {
                num_of_rows += dst_rate.size();
            }
        }
    }

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows),
                    event_names(num_of_rows), dst_epi_states(num_of_rows);
    NumericVector rates(num_of_rows), times(num_of_rows);

    size_t i{0};
    const auto &tissue = sim_ptr->tissue();
    for (const auto &[time, species_rate_updates] : rate_update_history) {
        for (const auto &[species_id, event_rate_updates] : species_rate_updates) {
            const auto &species = tissue.get_species(species_id);
            const std::string mutant_name = species.get_mutant_name();
            const auto epi_state = species.get_epistate_name();
            for (const auto &[event_name, dst_rate] : event_rate_updates) {
                for (const auto &[dst_id, rate] : dst_rate) {
                    times[i] = time;
                    mutant_names[i] = mutant_name.c_str();
                    epi_states[i] = encode_epigenetic_state(species);
                    event_names[i] = event_name.c_str();

                    if (get_event_id(event_name)==CellEventType::DEATH) {
                        dst_epi_states[i] = NA_STRING;
                    } else {
                        const auto &dst_species = tissue.get_species(dst_id);
                        dst_epi_states[i] = encode_epigenetic_state(dst_species);
                    }
                    rates[i] = rate;

                    ++i;
                }
            }
        }
    }

    if (sim_ptr->get_epigenetic_state_names().size()==0) {
        return DataFrame::create(_["time"] = times, _["mutant"] = mutant_names,
                                 _["event"] = event_names, _["rate"] = rates);
    }

    return DataFrame::create(_["time"] = times, _["mutant"] = mutant_names,
                             _["epistate"] = epi_states, _["event"] = event_names,
                             _["first child epistate"] = epi_states, _["rate"] = rates);
}

void TissueSimulation::place_cell(const std::string &species_name,
                                  const CLONES::Mutants::Evolutions::AxisPosition &x,
                                  const CLONES::Mutants::Evolutions::AxisPosition &y)
{
    if (sim_ptr->tissue().num_of_mutated_cells() > 0) {
        Rcpp::warning("Warning: the tissue already contains a cell.");
    }

    const auto &species = sim_ptr->tissue().get_species(species_name);

    sim_ptr->place_cell(species.get_id(), {x, y});
}

void TissueSimulation::place_cell(const SEXP &species_name,  const SEXP &x, const SEXP &y)
{
    const auto C_species_name = FromSEXP::get<std::string>(species_name, "1st parameter", "string");
    const auto C_x = FromSEXP::get<CLONES::Mutants::Evolutions::AxisPosition>(x, "2nd parameter", "natural number");
    const auto C_y = FromSEXP::get<CLONES::Mutants::Evolutions::AxisPosition>(y, "3rd parameter", "natural number");

    if (sim_ptr->tissue().num_of_mutated_cells() > 0) {
        Rcpp::warning("Warning: the tissue already contains a cell.");
    }

    const auto &species = sim_ptr->tissue().get_species(C_species_name);

    sim_ptr->place_cell(species.get_id(), {C_x, C_y});
}

Rcpp::List
TissueSimulation::get_cells(const CLONES::Mutants::Evolutions::Tissue &tissue) const
{
    namespace RS = CLONES::Mutants::Evolutions;

    std::vector<RS::AxisPosition> upper_corner = tissue.size();
    upper_corner.resize(2);

    for (auto &value : upper_corner) {
        --value;
    }

    return get_cells(tissue, {0, 0}, upper_corner);
}

Rcpp::List
TissueSimulation::get_cell(const CLONES::Mutants::Evolutions::Tissue &tissue,
                           const CLONES::Mutants::Evolutions::AxisPosition &x,
                           const CLONES::Mutants::Evolutions::AxisPosition &y) const
{
    const auto cell_proxy = tissue({x, y});

    return wrap_a_cell(cell_proxy);
}

Rcpp::List TissueSimulation::get_cells(
    const CLONES::Mutants::Evolutions::Tissue &tissue,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner) const
{
    std::set<CLONES::Mutants::SpeciesId> species_ids;

    for (const auto &species : tissue) {
        species_ids.insert(species.get_id());
    }

    return get_cells(tissue, lower_corner, upper_corner, species_ids,
                     tissue.get_epigenetic_state_names());
}

Rcpp::List TissueSimulation::get_cells(const SEXP &first_param,
                                       const SEXP &second_param) const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants::Evolutions;

    if (TYPEOF(first_param) != TYPEOF(second_param)) {
        std::ostringstream oss;

        oss << "The two parameters have different types: " << type2name(first_param)
            << " != " << type2name(second_param);

        Rcpp::stop(oss.str());
    }

    switch (TYPEOF(first_param)) {
    case INTSXP:
    case REALSXP:
    {
        return get_cells(as<std::vector<AxisPosition>>(first_param),
                         as<std::vector<AxisPosition>>(second_param));
    }
    case STRSXP:
    {
        return get_cells(as<std::vector<std::string>>(first_param),
                         as<std::vector<std::string>>(second_param));
    }
    default:
    {
        std::ostringstream oss;

        oss << "Invalid parameter type " << type2name(first_param);

        Rcpp::stop(oss.str());
    }
    }
}

Rcpp::List
TissueSimulation::get_cells(const std::vector<std::string> &species_filter,
                            const std::vector<std::string> &epigenetic_filter) const
{
    namespace RS = CLONES::Mutants::Evolutions;

    std::vector<RS::AxisPosition> upper_corner = sim_ptr->tissue().size();
    upper_corner.resize(2);

    for (auto &value : upper_corner) {
        --value;
    }

    return get_cells(sim_ptr->tissue(), {0, 0}, upper_corner, species_filter,
                     epigenetic_filter);
}

Rcpp::List TissueSimulation::get_cells(
    const CLONES::Mutants::Evolutions::Tissue &tissue,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner,
    const std::vector<std::string> &mutant_filter,
    const std::vector<std::string> &epigenetic_filter) const
{
    std::set<std::string> mutant_set(mutant_filter.begin(), mutant_filter.end());
    std::set<std::string> epigenetic_set(epigenetic_filter.begin(),
                                         epigenetic_filter.end());

    auto species_ids = get_species_ids_from_mutant_name(tissue, mutant_set);

    return get_cells(tissue, lower_corner, upper_corner, species_ids, epigenetic_set);
}

Rcpp::List TissueSimulation::get_counts() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    size_t num_of_rows = sim_ptr->tissue().num_of_species();

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
    IntegerVector counts(num_of_rows), overall(num_of_rows);

    size_t i{0};
    for (const auto &species : sim_ptr->tissue()) {
        mutant_names[i] = species.get_mutant_name();
        epi_states[i] = encode_epigenetic_state(species);
        counts[i] = species.num_of_cells();
        overall[i] = species.num_of_simulated_cells();
        ++i;
    }

    if (sim_ptr->tissue().get_epigenetic_state_names().size() == 0) {
        return DataFrame::create(_["mutant"] = mutant_names,
                                 _["counts"] = counts, _["overall"] = overall);
    }

    return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states,
                             _["counts"] = counts, _["overall"] = overall);
}

std::map<CLONES::Mutants::SpeciesId, std::string>
get_species_id2name(const CLONES::Mutants::Evolutions::Tissue &tissue)
{
    std::map<CLONES::Mutants::SpeciesId, std::string> id2name;
    for (const auto &species : tissue) {
        id2name[species.get_id()] = species.get_name();
    }

    return id2name;
}

Rcpp::List TissueSimulation::get_added_cells() const
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    namespace RS = CLONES::Mutants::Evolutions;

    size_t num_of_rows = sim_ptr->get_added_cells().size();

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
    IntegerVector position_x(num_of_rows), position_y(num_of_rows);
    NumericVector time(num_of_rows);

    size_t i{0};
    for (const auto &added_cell : sim_ptr->get_added_cells()) {
        const auto &species = sim_ptr->tissue().get_species(added_cell.species_id);
        mutant_names[i] = species.get_mutant_name();
        epi_states[i] = encode_epigenetic_state(species);
        position_x[i] = added_cell.x;
        position_y[i] = added_cell.y;
        time[i] = added_cell.time;
        ++i;
    }

    if (sim_ptr->tissue().get_epigenetic_state_names().size() == 0) {
        return DataFrame::create(_["mutant"] = mutant_names,
                                 _["position_x"] = position_x, _["position_y"] = position_y,
                                 _["time"] = time);
    }

    return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states,
                             _["position_x"] = position_x, _["position_y"] = position_y,
                             _["time"] = time);
}

// sorting LineageEdge by time
struct TimedLineageEdge : public CLONES::Mutants::Evolutions::LineageEdge
{
    CLONES::Time time;

    TimedLineageEdge() : CLONES::Mutants::Evolutions::LineageEdge(), time(0) {}

    TimedLineageEdge(const CLONES::Mutants::Evolutions::LineageEdge &edge,
                     const CLONES::Time &time)
        : CLONES::Mutants::Evolutions::LineageEdge(edge), time(time)
    {}
};

struct TimedLineageEdgeCmp
{
    bool operator()(const TimedLineageEdge &a, const TimedLineageEdge &b)
    {
        return (a.time < b.time ||
                (a.time == b.time && (a.get_ancestor() < b.get_ancestor())) ||
                (a.time == b.time && (a.get_ancestor() == b.get_ancestor()) &&
                 (a.get_progeny() < b.get_progeny())));
    }
};

std::vector<TimedLineageEdge>
sorted_timed_edges(const CLONES::Mutants::Evolutions::TissueSimulation &simulation)
{
    const auto &lineage_graph = simulation.get_lineage_graph();
    const size_t num_of_edges = lineage_graph.num_of_edges();

    std::vector<TimedLineageEdge> timed_edges;

    timed_edges.reserve(num_of_edges);

    for (const auto &[edge, edge_time] : lineage_graph) {
        timed_edges.push_back({edge, edge_time});
    }

    TimedLineageEdgeCmp cmp;
    sort(timed_edges.begin(), timed_edges.end(), cmp);

    return timed_edges;
}

Rcpp::List TissueSimulation::get_lineage_graph() const
{
    using namespace Rcpp;
    const auto species_id2name = get_species_id2name(sim_ptr->tissue());

    const auto timed_edges = sorted_timed_edges(*sim_ptr);

    CharacterVector ancestors(timed_edges.size()), progeny(timed_edges.size());
    NumericVector first_cross(timed_edges.size());

    size_t i{0};
    for (const auto &timed_edge : timed_edges) {
        ancestors[i] = (timed_edge.get_ancestor() != WILD_TYPE_SPECIES
                            ? species_id2name.at(timed_edge.get_ancestor())
                            : "Wild-type");

        progeny[i] = (timed_edge.get_progeny() != WILD_TYPE_SPECIES
                          ? species_id2name.at(timed_edge.get_progeny())
                          : "Wild-type");
        first_cross[i] = timed_edge.time;

        ++i;
    }

    return DataFrame::create(_["ancestor"] = ancestors, _["progeny"] = progeny,
                             _["first_cross"] = first_cross);
}

inline void validate_non_empty_tissue(const CLONES::Mutants::Evolutions::Tissue &tissue)
{
    if (tissue.num_of_cells() == 0) {
        Rcpp::stop("The tissue does not contain any cell.");
    }
}

void TissueSimulation::run_up_to_time(const CLONES::Time &time)
{
    run_up_to_time(time, false);
}

void TissueSimulation::run_up_to_time(const CLONES::Time &time, const bool quiet)
{
    validate_non_empty_tissue(sim_ptr->tissue());

    RTest<CLONES::Mutants::Evolutions::TimeTest> ending_test{time};

    CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    sim_ptr->run(ending_test, progress_bar);
}

void TissueSimulation::run_up_to_size(const std::string &species_name,
                                      const size_t &num_of_cells)
{
    run_up_to_size(species_name, num_of_cells, false);
}

void TissueSimulation::run_up_to_size(const std::string &species_name,
                                      const size_t &num_of_cells, const bool quiet)
{
    validate_non_empty_tissue(sim_ptr->tissue());

    const auto &species_id = sim_ptr->tissue().get_species(species_name).get_id();

    RTest<CLONES::Mutants::Evolutions::SpeciesCountTest> ending_test{species_id,
                                                                    num_of_cells};

    CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    sim_ptr->run(ending_test, progress_bar);
}

void TissueSimulation::run_up_to_event(const std::string &event,
                                       const std::string &species_name,
                                       const size_t &num_of_events)
{
    run_up_to_event(event, species_name, num_of_events, false);
}

void TissueSimulation::run_up_to_event(const std::string &event,
                                       const std::string &species_name,
                                       const size_t &num_of_events, const bool quiet)
{
    validate_non_empty_tissue(sim_ptr->tissue());

    namespace RS = CLONES::Mutants::Evolutions;

    const auto &species_id = sim_ptr->tissue().get_species(species_name).get_id();

    RTest<RS::EventCountTest> ending_test{get_event_id(event), species_id,
                                          num_of_events};

    CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    sim_ptr->run(ending_test, progress_bar);
}

void TissueSimulation::run_until(const Logics::Formula &formula)
{
    run_until(formula, false);
}

void TissueSimulation::run_until(const Logics::Formula &formula, const bool quiet)
{
    validate_non_empty_tissue(sim_ptr->tissue());

    RTest<CLONES::Mutants::Evolutions::FormulaTest> ending_test{formula};

    CLONES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

    sim_ptr->run(ending_test, progress_bar);
}

Rcpp::List TissueSimulation::get_firings() const
{
    using namespace Rcpp;

    const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

    auto df = get_firing_history(last_time_sample, last_time_sample);

    if (sim_ptr->tissue().get_epigenetic_state_names().size() == 0) {
        return DataFrame::create(_["event"] = df["event"], _["mutant"] = df["mutant"],
                                 _["fired"] = df["fired"]);
    }

    return DataFrame::create(_["event"] = df["event"], _["mutant"] = df["mutant"],
                             _["epistate"] = df["epistate"], _["fired"] = df["fired"]);
}

Rcpp::List TissueSimulation::get_firing_history(const CLONES::Time &minimum_time) const
{
    if (sim_ptr->get_statistics().get_history().size() == 0) {
        return get_firing_history(0, 0);
    }

    const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

    return get_firing_history(minimum_time, last_time_sample);
}

size_t TissueSimulation::count_history_sample_in(const CLONES::Time &minimum_time,
                                                 const CLONES::Time &maximum_time) const
{
    size_t num_of_samples{0};
    const auto &history = sim_ptr->get_statistics().get_history();
    auto series_it = history.lower_bound(minimum_time);
    while (series_it != history.end() && series_it->first <= maximum_time) {
        ++num_of_samples;
        ++series_it;
    }

    return num_of_samples;
}

Rcpp::List TissueSimulation::get_firing_history(const CLONES::Time &minimum_time,
                                                const CLONES::Time &maximum_time) const
{
    using namespace Rcpp;

    const bool with_epigenetics{sim_ptr->tissue().get_epigenetic_state_names().size()>0};

    const size_t num_of_events = CLONES::Mutants::cell_event_names.size();
    const size_t rows_per_sample = (num_of_events-(with_epigenetics?1:2)) * sim_ptr->tissue().num_of_species();
    const size_t num_of_rows =
        count_history_sample_in(minimum_time, maximum_time) * rows_per_sample;

    CharacterVector events(num_of_rows), mutant_names(num_of_rows),
        epi_states(num_of_rows);
    IntegerVector firings(num_of_rows);
    NumericVector times(num_of_rows);

    size_t i{0};
    const auto &history = sim_ptr->get_statistics().get_history();
    auto series_it = history.lower_bound(minimum_time);
    while (series_it != history.end() && series_it->first <= maximum_time) {
        const auto &time = series_it->first;
        const auto &t_stats = series_it->second;
        for (const auto &species : sim_ptr->tissue()) {
            for (const auto &[event_code, event_name] : CLONES::Mutants::cell_event_names) {
                if (event_code != CLONES::Mutants::CellEventType::MUTATION
                        && (with_epigenetics
                            || event_code != CLONES::Mutants::CellEventType::DUP_AND_EPI_SWITCH)) {
                    events[i] = event_name;
                    mutant_names[i] = species.get_mutant_name();
                    epi_states[i] = encode_epigenetic_state(species);

                    const auto &species_it = t_stats.find(species.get_id());
                    if (species_it != t_stats.end()) {
                        firings[i] = count_events(species_it->second, event_code);
                    } else {
                        firings[i] = 0;
                    }
                    times[i] = time;
                    ++i;
                }
            }
        }
        ++series_it;
    }

    if (sim_ptr->tissue().get_epigenetic_state_names().size() == 0) {
        return DataFrame::create(_["event"] = events, _["mutant"] = mutant_names,
                                 _["fired"] = firings, _["time"] = times);
    }

    return DataFrame::create(_["event"] = events, _["mutant"] = mutant_names,
                             _["epistate"] = epi_states, _["fired"] = firings,
                             _["time"] = times);
}

Rcpp::List TissueSimulation::get_count_history(const CLONES::Time &minimum_time) const
{
    if (sim_ptr->get_statistics().get_history().size() == 0) {
        return get_count_history(0, 0);
    }

    const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

    return get_count_history(minimum_time, last_time_sample);
}

Rcpp::List TissueSimulation::get_count_history(const CLONES::Time &minimum_time,
                                               const CLONES::Time &maximum_time) const
{
    using namespace Rcpp;

    const size_t rows_per_sample = sim_ptr->tissue().num_of_species();
    const size_t num_of_rows =
        count_history_sample_in(minimum_time, maximum_time) * rows_per_sample;

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
    IntegerVector counts(num_of_rows);
    NumericVector times(num_of_rows);

    bool with_epigenetics{false};

    size_t i{0};
    const auto &history = sim_ptr->get_statistics().get_history();
    auto series_it = history.lower_bound(minimum_time);
    while (series_it != history.end() && series_it->first <= maximum_time) {
        const auto &time = series_it->first;
        const auto &t_stats = series_it->second;
        for (const auto &species : sim_ptr->tissue()) {
            mutant_names[i] = species.get_mutant_name();

            const auto epistate = species.get_epistate_name();
            if (epistate != "") {
                epi_states[i] = epistate;
                with_epigenetics = true;
            }

            const auto &species_it = t_stats.find(species.get_id());
            if (species_it != t_stats.end()) {
                counts[i] = species_it->second.curr_cells;
            } else {
                counts[i] = 0;
            }
            times[i] = time;
            ++i;
        }
        ++series_it;
    }

    if (with_epigenetics) {
        return DataFrame::create(_["mutant"] = mutant_names, _["epistate"] = epi_states,
                                _["count"] = counts, _["time"] = times);
    }
    return DataFrame::create(_["mutant"] = mutant_names, _["count"] = counts,
                             _["time"] = times);
}

Rcpp::IntegerVector TissueSimulation::get_tissue_size() const
{
    auto size_vect = sim_ptr->tissue().size();

    return {size_vect[0], size_vect[1]};
}

void TissueSimulation::set_rate(const std::string &species_name,
                                const std::string &event_name,
                                const double &rate)
{
    using namespace CLONES::Mutants;

    if ((event_name != get_event_name(CellEventType::DEATH))
            && (event_name != get_event_name(CellEventType::DUPLICATION))) {

        Rcpp::stop(("set_rate(species_name, event_name, rate) can only be used for "
                    "setting the rates of ")
                   + get_event_name(CellEventType::DEATH)
                   + " and "
                   + get_event_name(CellEventType::DUPLICATION));
    }

    set_rate(species_name, event_name, species_name, rate);
}

void TissueSimulation::set_rate(const std::string &src_species_name,
                                const std::string &event_name,
                                std::string dst_name,
                                const double &rate)
{
    try {
        auto &src_species = sim_ptr->tissue().get_species(src_species_name);

        if (sim_ptr->tissue().knowns_epigenetic_state(dst_name)) {
            using namespace CLONES::Mutants;

            dst_name = SpeciesName{src_species.get_mutant_name(), dst_name};
        } else {
            if (!sim_ptr->tissue().knowns_species(dst_name)) {
                Rcpp::stop(("The third parameter must be either a species "
                            "or an epigenetic state. \"") + dst_name
                            + "\" is none of the two.");
            }
        }

        auto &dst_species = sim_ptr->tissue().get_species(dst_name);

        src_species.set_rate(get_event_id(event_name), dst_species, rate);

        auto &species_update = rate_update_history[sim_ptr->get_time()][src_species.get_id()];

        species_update[event_name][dst_species.get_id()] = rate;
    } catch (const std::exception& ex) {
        Rcpp::stop(ex.what());
    }

    // save changes
    CLONES::Archive::Binary::Out ruh_archive(get_rates_update_history_path());

    ruh_archive & rate_update_history;
}

CLONES::Mutants::SpeciesName TissueSimulation::get_species_name(const std::string& mutant_name, const std::string& species_name) const
{
    using namespace CLONES::Mutants;

    if (!sim_ptr->tissue().knowns_species(species_name)) {
        if (!sim_ptr->tissue().knowns_epigenetic_state(species_name)) {
            Rcpp::stop("\"" + species_name + "\" is neither a species nor an epigenetic state.");
        }

        return SpeciesName{mutant_name, species_name};
    }

    SpeciesName sname{species_name};

    if (sname.get_mutant_name() != mutant_name) {
        Rcpp::stop("\"" + species_name + "\" has a mutant different from \"" + mutant_name + "\"");
    }

    return sname;
}

template<typename T>
Rcpp::List extract_list(T value, const std::string& what)
{
    using namespace Rcpp;

    RObject obj{value};

    if (!Rf_isNewList(obj)) {
        Rcpp::stop(what + " must is a named list.");
    }

    List list = as<List>(obj);

    if (!list.hasAttribute("names")) {
        Rcpp::stop(what + " must is a named list.");
    }

    return list;
}

void TissueSimulation::set_species_rates(const std::string &species_name,
                                         const Rcpp::List &rates)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    if (!rates.hasAttribute("names")) {
        Rcpp::stop("The species rate specification is supposed to be a named list.");
    }

    SpeciesName sname{species_name};

    CharacterVector nv = rates.names();
    for (int i = 0; i < nv.size(); ++i){
        const std::string event_name = as<std::string>(nv[i]);
        if (event_name == get_event_name(CellEventType::DEATH)
                || event_name == get_event_name(CellEventType::DUPLICATION)) {

            const double rate_value = rates[i];
            set_rate(species_name, event_name, rate_value);
        } else { // the event must be a switch and the name may be either an
                 // epigenetic state or a species

            auto dst_species_name = get_species_name(sname.get_mutant_name(), event_name);
            const double rate_value = rates[i];

            set_rate(species_name, get_event_name(CellEventType::DUP_AND_EPI_SWITCH),
                     dst_species_name, rate_value);
        }
    }
}

std::set<std::string> collect_string_set(const Rcpp::CharacterVector &vector)
{
    std::set<std::string> strings;

    for (auto i=0; i< vector.size(); ++i) {
        strings.insert(as<std::string>(vector[i]));
    }

    return strings;
}

bool has_df_epi_columns(const std::set<std::string>& names,
                        const std::list<std::string> epi_columns)
{
    std::set<std::string> missing;
    for (const std::string& col_name : epi_columns) {
        if (names.count(col_name)==0) {
            missing.insert(col_name);
        }
    }

    if (missing.size()!=0 && missing.size()!=epi_columns.size()) {
        Rcpp::stop(("The data frame must either contains both the columns "
                    "\"epistate\" or \"first.child.epistate\" or none of them. "
                    "It misses the column \"") + *(missing.begin()) + "\".");
    }

    return missing.size()==epi_columns.size();
}

void validate_column_names(const Rcpp::DataFrame &df, bool has_epistate)
{
    using namespace Rcpp;

    CharacterVector R_col_names = df.attr("names");

    std::set<std::string> df_names;
    for (auto i=0; i< R_col_names.size(); ++i) {
        df_names.emplace(R_col_names[i]);
    }

    std::list<std::string> columns{"mutant", "event", "rate"};
    for (const auto& col_name : columns) {
        if (df_names.count(col_name)==0) {
            Rcpp::stop("The data fra me misses the column \"" + col_name + "\".");
        }
    }

    const std::list<std::string> epi_columns{"epistate", "first.child.epistate"};

    for (const auto& name : epi_columns) {
        if (has_epistate) {
            if (df_names.count(name)==0) {
                Rcpp::stop(("The simulation contains epigenetic states and the data"
                            " frame misses the columns \"") + name + "\".");
            }
        } else {
            if (df_names.count(name)>0) {
                Rcpp::stop(("The simulation has no epigenetic state, while the data"
                            " frame contains the column \"") + name + "\".");
            }
        }
    }
}


void TissueSimulation::set_rate_from_df_row(const SEXP& mutant,
                                            const SEXP& event_name,
                                            const double& rate, const size_t& row_num)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    try {
        set_rate(as<std::string>(mutant), as<std::string>(event_name), rate);
    } catch(const std::exception &ex) {
        stop("DataFrame line " + std::to_string(row_num) + ": " + ex.what());
    }
}

void TissueSimulation::set_rate_from_df_row(const SEXP& mutant,
                                            const SEXP& epistate,
                                            const SEXP& event_name,
                                            const SEXP& fc_epistate,
                                            const double& rate, const size_t& row_num)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    try {
        if (CharacterVector::is_na(epistate)) {
            set_rate(as<std::string>(mutant), as<std::string>(event_name), rate);
        } else {
            const std::string src_epistate = as<std::string>(epistate);

            std::string species = SpeciesName{as<std::string>(mutant), src_epistate};

            if (CharacterVector::is_na(fc_epistate)) {
                set_rate(species, as<std::string>(event_name), rate);
            } else {
                set_rate(species, as<std::string>(event_name), as<std::string>(fc_epistate), rate);
            }
        }
    } catch(const std::exception &ex) {
        stop("DataFrame line " + std::to_string(row_num) + ": " + ex.what());
    }
}

void TissueSimulation::set_rates(const Rcpp::DataFrame &rates)
{
    using namespace Rcpp;

    bool has_epistates = sim_ptr->tissue().get_epigenetic_state_names().size()>0;

    validate_column_names(rates, has_epistates);

    CharacterVector mutant_col = rates["mutant"];
    CharacterVector event_col = rates["event"];
    NumericVector rate_col = rates["rate"];

    if (has_epistates) {
        CharacterVector epistate_col = rates["epistate"];
        CharacterVector fc_epistate_col = rates["first.child.epistate"];

        for (auto i=0; i< rates.nrows(); ++i) {
            set_rate_from_df_row(mutant_col[i], epistate_col[i], event_col[i],
                                 fc_epistate_col[i], rate_col[i], i);
        }
    } else {
        for (auto i=0; i< rates.nrows(); ++i) {
            set_rate_from_df_row(mutant_col[i], event_col[i], rate_col[i], i);
        }
    }
}

void TissueSimulation::set_rates(const Rcpp::List &rates)
{
    using namespace Rcpp;

    if (is<DataFrame>(rates)) {
        DataFrame df = as<DataFrame>(rates);

        set_rates(df);

        return;
    }

    if (!rates.hasAttribute("names")) {
        Rcpp::stop("The parameter must be a named list whose names are species.");
    }

    CharacterVector nv = rates.names();
    for (int i = 0; i < nv.size(); ++i){
        const std::string species_name = as<std::string>(nv[i]);
        if (!sim_ptr->tissue().knowns_species(species_name)) {
            Rcpp::stop("\"" + species_name + "\" is not a species name.");
        }

        RObject obj{rates[i]};

        if (!Rf_isNewList(obj)) {
            Rcpp::stop("The rate specifications for \"" + species_name
                       + "\" must is a named list.");
        }

        CLONES::Mutants::SpeciesName sname{species_name};
        if (sname.get_epistate_name()=="") {
            set_rates(sname.get_mutant_name(), rates[i]);
        } else {
            List rate_list = List::create(
                _[sname.get_epistate_name()] = obj
            );

            set_rates(sname.get_mutant_name(), rate_list);
        }
    }
}

void TissueSimulation::set_rates(const std::string &mutant_name,
                                 const Rcpp::List &rates)
{
    using namespace Rcpp;
    using namespace CLONES::Mutants;

    if (!sim_ptr->knowns_mutant(mutant_name)) {
        if (!sim_ptr->knowns_species(mutant_name)) {
            Rcpp::stop("\"" + mutant_name
                       + "\" is neither a mutant nor a species.");
        }

        set_species_rates(mutant_name, rates);
    }

    auto m_view = sim_ptr->tissue().get_mutant_view(mutant_name);

    if (!rates.hasAttribute("names")) {
        Rcpp::stop("The second parameter must be a Rcpp::List "
                   "with the names attribute");
    }

    CharacterVector nv = rates.names();
    for (int i = 0; i < nv.size(); ++i){
        const std::string event_name = as<std::string>(nv[i]);
        if (event_name == get_event_name(CellEventType::DEATH)
                || event_name == get_event_name(CellEventType::DUPLICATION)) {
            for (auto& species: m_view) {
                const double rate_value = rates[i];
                set_rate(species.get_name(), event_name, rate_value);
            }
        }
    }

    for (int i = 0; i < nv.size(); ++i){
        std::string species_name = as<std::string>(nv[i]);

        if (species_name != get_event_name(CellEventType::DEATH)
                && species_name != get_event_name(CellEventType::DUPLICATION)) {
            species_name = get_species_name(mutant_name, species_name);

            auto epistate_rates = extract_list(rates[i], "The species rate specifications");

            set_species_rates(species_name, epistate_rates);
        }
    }
}

std::list<std::string> collect_names(const SEXP &names)
{
    using namespace Rcpp;

    if (is<String>(names)) {
        return {as<std::string>(names)};
    }

    if (is<List>(names) || is<CharacterVector>(names)) {
        return collect_strings(as<List>(names));
    }

    stop("The first parameter of must be either a mutant/species "
         "name or a list of mutant/species names.");
}

Rcpp::List TissueSimulation::choose_cell_in(const SEXP &names)
{
    const auto C_names = collect_names(names);

    const auto &cell = sim_ptr->choose_cell_in(C_names);

    return wrap_a_cell(cell);
}

Rcpp::List TissueSimulation::choose_cell_in(const SEXP &names,
        const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
        const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner)
{
    const auto rectangle = get_rectangle(lower_corner, upper_corner);

    const auto C_names = collect_names(names);

    const auto &cell = sim_ptr->choose_cell_in(C_names, rectangle);

    return wrap_a_cell(cell);
}

Rcpp::List TissueSimulation::choose_border_cell_in(const SEXP &names)
{
    const auto C_names = collect_names(names);

    const auto &cell = sim_ptr->choose_border_cell_in(C_names);

    return wrap_a_cell(cell);
}

Rcpp::List TissueSimulation::choose_border_cell_in(const SEXP &names,
        const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
        const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner)
{
    const auto rectangle = get_rectangle(lower_corner, upper_corner);

    const auto C_names = collect_names(names);

    const auto &cell = sim_ptr->choose_border_cell_in(C_names, rectangle);

    return wrap_a_cell(cell);
}

void TissueSimulation::mutate_progeny(const CLONES::Mutants::Evolutions::AxisPosition &x,
                                      const CLONES::Mutants::Evolutions::AxisPosition &y,
                                      const std::string &mutated_mutant)
{
    auto pos_in_tissue = get_position_in_tissue({x, y});

    namespace RS = CLONES::Mutants::Evolutions;

    sim_ptr->simulate_mutation(pos_in_tissue, mutated_mutant);
}

void TissueSimulation::mutate_progeny(const Rcpp::List &cell_position,
                                      const std::string &mutated_mutant)
{
    using namespace Rcpp;

    namespace RS = CLONES::Mutants::Evolutions;

    std::vector<CLONES::Mutants::Evolutions::AxisPosition> vector_position;

    for (const std::string axis : {"x", "y"}) {
        auto field = "position_" + axis;
        if (!cell_position.containsElementNamed(field.c_str())) {
            Rcpp::stop("Missing \"" + field + "\" element from the Rcpp::List.");
        }
        vector_position.push_back(as<RS::AxisPosition>(cell_position[field]));
    }

    return mutate_progeny(vector_position[0], vector_position[1], mutated_mutant);
}

bool TissueSimulation::already_collected_sample(const std::string &sample_name) const
{
    for (const auto &sample : sim_ptr->get_tissue_samples()) {
        if (sample.get_name() == sample_name) {
            return true;
        }
    }

    return false;
}

void TissueSimulation::validate_usable_sample_name(const std::string &sample_name) const
{
    if (sample_name == "normal sample" || sample_name == "normal.sample") {
        Rcpp::stop("Sample name \"normal sample\" is reserved.");
    }

    if (already_collected_sample(sample_name)) {
        Rcpp::stop("The sample \"" + sample_name + "\" has been already collected.");
    }
}

void TissueSimulation::sample_cells(
    const std::string &sample_name,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner) const
{
    using namespace CLONES::Mutants;

    const auto bounding_box = get_rectangle(lower_corner, upper_corner);
    const auto num_of_cells = bounding_box.size();

    sample_cells(sample_name, lower_corner, upper_corner, num_of_cells);
}

void TissueSimulation::sample_cells(const std::string &sample_name,
                                    const size_t &num_of_cells) const
{
    std::vector<CLONES::Mutants::Evolutions::AxisPosition> lower_corner, upper_corner;

    for (const auto axis_size : sim_ptr->tissue().size()) {
        lower_corner.push_back(0);
        upper_corner.push_back(axis_size - 1);
    }

    sample_cells(sample_name, lower_corner, upper_corner, num_of_cells);
}

void TissueSimulation::sample_cells(
    const std::string &sample_name,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &lower_corner,
    const std::vector<CLONES::Mutants::Evolutions::AxisPosition> &upper_corner,
    const size_t &num_of_cells) const
{
    using namespace CLONES::Mutants;

    validate_usable_sample_name(sample_name);
    save_tissue(pre_sample_tissue_image_path(sample_name));

    auto bounding_box = get_rectangle(lower_corner, upper_corner);

    CLONES::Mutants::Evolutions::SampleSpecification spec(sample_name, bounding_box,
                                                         num_of_cells);

    sim_ptr->sample_tissue(spec);
}

std::map<CLONES::Mutants::SpeciesId, size_t>
count_cells_in(const CLONES::Mutants::Evolutions::Tissue &tissue, const uint16_t &init_x,
               const uint16_t &init_y, const uint16_t &width, const uint16_t &height)
{
    std::map<CLONES::Mutants::SpeciesId, size_t> counter;
    auto sizes = tissue.size();

    uint16_t x_max = std::min(static_cast<uint16_t>(init_x + width), sizes[0]);
    uint16_t y_max = std::min(static_cast<uint16_t>(init_y + height), sizes[1]);

    for (uint16_t x = init_x; x < x_max; ++x) {
        for (uint16_t y = init_y; y < y_max; ++y) {
            auto cell_proxy = tissue({x, y});
            if (!cell_proxy.is_wild_type()) {
                const CLONES::Mutants::Evolutions::CellInTissue &cell = cell_proxy;

                const auto species_id = cell.get_species_id();

                auto found = counter.find(species_id);

                if (found == counter.end()) {
                    counter.insert({species_id, 1});
                } else {
                    ++(found->second);
                }
            }
        }
    }

    return counter;
}

inline std::map<CLONES::Mutants::SpeciesId, size_t>
count_cells_in(const CLONES::Mutants::Evolutions::Tissue &tissue,
               const TissueRectangle &tumour_bounding_box, const uint16_t &grid_x,
               const uint16_t &grid_y, const uint16_t &width, const uint16_t &height)
{
    const uint16_t x = grid_x * width + tumour_bounding_box.lower_corner.x,
                   y = grid_y * height + tumour_bounding_box.lower_corner.y;

    return count_cells_in(tissue, x, y, width, height);
}

inline TissueRectangle get_tissue_rectangle(const TissueRectangle &tumour_bounding_box,
                                            const uint16_t &grid_x,
                                            const uint16_t &grid_y, const uint16_t &width,
                                            const uint16_t &height)
{
    using namespace CLONES::Mutants::Evolutions;

    const uint16_t x = grid_x * width + tumour_bounding_box.lower_corner.x,
                   y = grid_y * height + tumour_bounding_box.lower_corner.y;

    return TissueRectangle(PositionInTissue{x, y}, width, height);
}

std::set<CLONES::Mutants::SpeciesId>
collect_species_of(const CLONES::Mutants::Evolutions::TissueSimulation &simulation,
                   const std::string &mutant_name)
{
    std::set<CLONES::Mutants::SpeciesId> species_ids;

    if (mutant_name.back() != '+' && mutant_name.back() != '-') {
        const auto &tissue = simulation.tissue();

        auto mutant_id = simulation.find_mutant_id(mutant_name);

        for (const auto &species : tissue) {
            if (species.get_mutant_id() == mutant_id) {
                species_ids.insert(species.get_id());
            }
        }
    }

    return species_ids;
}

TissueRectangle TissueSimulation::get_tumour_bounding_box() const
{
    using namespace CLONES::Mutants::Evolutions;
    const auto &tissue = sim_ptr->tissue();
    const auto tissue_sizes = tissue.size();

    PositionInTissue lower_corner(static_cast<AxisSize>(tissue_sizes[0]),
                                  static_cast<AxisSize>(tissue_sizes[1])),
        upper_corner{0, 0};

    for (uint16_t grid_x = 0; grid_x < tissue_sizes[0]; ++grid_x) {
        for (uint16_t grid_y = 0; grid_y < tissue_sizes[1]; ++grid_y) {
            CLONES::Mutants::Evolutions::PositionInTissue pos{grid_x, grid_y};
            if (!tissue(pos).is_wild_type()) {
                if (grid_x < lower_corner.x) {
                    lower_corner.x = grid_x;
                }
                if (grid_y < lower_corner.y) {
                    lower_corner.y = grid_y;
                }
                if (grid_x > upper_corner.x) {
                    upper_corner.x = grid_x;
                }
                if (grid_y > upper_corner.y) {
                    upper_corner.y = grid_y;
                }
            }
        }
    }

    return {lower_corner, upper_corner};
}

template <typename T> inline T div_ceil(const T &dividend, const T &divisor)
{
    if (dividend == 0)
        return dividend;

    return 1 + (dividend - 1) / divisor;
}

struct SpeciesConstraint
{
    const std::string species_name;
    const CLONES::Mutants::SpeciesId species_id;
    const size_t min_num_of_cells;

    SpeciesConstraint(const std::string &species_name,
                      const CLONES::Mutants::SpeciesId &species_id,
                      const size_t min_num_of_cells)
        : species_name(species_name), species_id(species_id),
          min_num_of_cells(min_num_of_cells)
    {}

    bool
    is_satified(const std::map<CLONES::Mutants::SpeciesId, size_t> &num_of_cells) const
    {
        auto found = num_of_cells.find(species_id);

        if (found != num_of_cells.end()) {
            return found->second >= min_num_of_cells;
        }

        return false;
    }
};

std::list<SpeciesConstraint>
get_species_constraints(const CLONES::Mutants::Evolutions::TissueSimulation &simulation,
                        const Rcpp::IntegerVector &minimum_cell_vector)
{
    std::list<SpeciesConstraint> species_constraints;

    const auto &tissue = simulation.tissue();

    const Rcpp::CharacterVector names = minimum_cell_vector.names();
    for (auto i = 0; i < minimum_cell_vector.size(); ++i) {
        const auto name = Rcpp::as<std::string>(names[i]);

        if (name.back() == '+' || name.back() == '-') {
            bool name_found{false};
            for (const auto &species : tissue) {
                if (species.get_name() == name) {
                    name_found = true;
                    const auto &threshold = minimum_cell_vector[i];

                    if (minimum_cell_vector[i] < 0) {
                        Rcpp::stop("The minimum number of cells must be "
                                   "a non-negative number. Specified " +
                                   std::to_string(threshold) + " for species \"" + name +
                                   "\".");
                    }
                    species_constraints.emplace_back(name, species.get_id(),
                                                     static_cast<size_t>(threshold));
                }
            }

            if (!name_found) {
                throw std::out_of_range("Unknown species \"" + name + "\"");
            }
        } else {
            simulation.find_mutant_id(name);
        }
    }

    return species_constraints;
}

struct MutantConstraint
{
    const std::string mutant_name;
    const std::set<CLONES::Mutants::SpeciesId> species_ids;
    const size_t min_num_of_cells;

    MutantConstraint(const std::string &mutant_name,
                     const std::set<CLONES::Mutants::SpeciesId> &species_ids,
                     const size_t &min_num_of_cells)
        : mutant_name(mutant_name), species_ids(species_ids),
          min_num_of_cells(min_num_of_cells)
    {}

    MutantConstraint(const std::string &mutant_name,
                     std::set<CLONES::Mutants::SpeciesId> &&species_ids,
                     size_t &&min_num_of_cells)
        : mutant_name(mutant_name), species_ids(std::move(species_ids)),
          min_num_of_cells(min_num_of_cells)
    {}

    bool
    is_satified(const std::map<CLONES::Mutants::SpeciesId, size_t> &num_of_cells) const
    {
        size_t total{0};

        for (const auto &species_id : species_ids) {
            auto found = num_of_cells.find(species_id);

            if (found != num_of_cells.end()) {
                total += found->second;
            }
        }

        return total > min_num_of_cells;
    }
};

std::list<MutantConstraint>
get_mutant_constraints(const CLONES::Mutants::Evolutions::TissueSimulation &simulation,
                       const Rcpp::IntegerVector &minimum_cell_vector)
{
    std::list<MutantConstraint> mutant_constraints;

    const Rcpp::CharacterVector names = minimum_cell_vector.names();
    for (auto i = 0; i < minimum_cell_vector.size(); ++i) {
        const std::string name = Rcpp::as<std::string>(names[i]);

        auto species_ids = collect_species_of(simulation, name);

        if (species_ids.size() > 0) {
            const auto &threshold = minimum_cell_vector[i];

            if (threshold < 0) {
                Rcpp::stop("The minimum number of cells must be "
                           "a non-negative number. Specified " +
                           std::to_string(threshold) + " for species \"" + name + "\".");
            }
            mutant_constraints.emplace_back(name, std::move(species_ids),
                                            static_cast<size_t>(threshold));
        }
    }

    return mutant_constraints;
}

bool constraints_satisfied(const std::map<CLONES::Mutants::SpeciesId, size_t> &cell_counts,
                           const std::list<SpeciesConstraint> &species_constraints,
                           const std::list<MutantConstraint> &mutant_constraints)
{
    for (const auto &constraint : species_constraints) {
        if (!constraint.is_satified(cell_counts)) {
            return false;
        }
    }

    for (const auto &constraint : mutant_constraints) {
        if (!constraint.is_satified(cell_counts)) {
            return false;
        }
    }

    return true;
}

std::vector<TissueRectangle>
TissueSimulation::find_all_samples(const Rcpp::IntegerVector &minimum_cell_vector,
                                   const uint16_t &width, const uint16_t &height) const
{
    auto species_constraints = get_species_constraints(*sim_ptr, minimum_cell_vector);
    auto mutant_constraints = get_mutant_constraints(*sim_ptr, minimum_cell_vector);

    auto t_bbox = get_tumour_bounding_box();

    const auto &tissue = sim_ptr->tissue();
    const auto t_width = t_bbox.upper_corner.x - t_bbox.lower_corner.x;
    const auto t_height = t_bbox.upper_corner.y - t_bbox.lower_corner.y;

    const uint16_t grid_width = t_width / width + ((t_width % width > 0 ? 1 : 0));
    const uint16_t grid_height = t_height / height + ((t_height % height > 0 ? 1 : 0));

    std::vector<TissueRectangle> rectangles;
    for (size_t grid_x = 0; grid_x < grid_width; ++grid_x) {
        for (size_t grid_y = 0; grid_y < grid_height; ++grid_y) {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);

            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                rectangles.push_back(
                    get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height));
            }
        }
    }

    return rectangles;
}

std::vector<TissueRectangle>
TissueSimulation::search_samples(const Rcpp::IntegerVector &minimum_cell_vector,
                                 const uint16_t &width, const uint16_t &height,
                                 const size_t num_of_samples, const int seed) const
{
    auto rectangles = find_all_samples(minimum_cell_vector, width, height);
    std::mt19937 gen(seed);

    if (rectangles.size() < num_of_samples) {
        std::ostringstream oss;
        if (rectangles.size() == 0) {
            oss << "No box satisfies";
        } else if (rectangles.size() == 1) {
            oss << "Only 1 sample satisfies";
        } else {
            oss << "Only " << rectangles.size() << " samples satisfy";
        }
        Rcpp::stop(oss.str() + " the constraints.");
    }

    std::vector<TissueRectangle> output;
    while (output.size() < num_of_samples) {
        std::uniform_int_distribution<size_t> selector(0, rectangles.size() - 1);

        size_t pos = selector(gen);

        output.push_back(rectangles[pos]);
        std::swap(rectangles[pos], rectangles[rectangles.size() - 1]);
        rectangles.resize(rectangles.size() - 1);
    }

    return output;
}

TissueRectangle
TissueSimulation::search_sample(const Rcpp::IntegerVector &minimum_cell_vector,
                                const uint16_t &width, const uint16_t &height) const
{
    auto species_constraints = get_species_constraints(*sim_ptr, minimum_cell_vector);
    auto mutant_constraints = get_mutant_constraints(*sim_ptr, minimum_cell_vector);

    auto t_bbox = get_tumour_bounding_box();

    const auto &tissue = sim_ptr->tissue();
    const auto t_width = t_bbox.upper_corner.x - t_bbox.lower_corner.x;
    const auto t_height = t_bbox.upper_corner.y - t_bbox.lower_corner.y;

    const uint16_t grid_width = t_width / width + ((t_width % width > 0 ? 1 : 0));
    const uint16_t grid_height = t_height / height + ((t_height % height > 0 ? 1 : 0));

    const uint16_t diag_size =
        div_ceil(std::min(grid_width, grid_height), static_cast<uint16_t>(2));

    for (uint16_t diag = 0; diag < diag_size; ++diag) {
        uint16_t grid_x = diag, grid_y = diag;

        for (; grid_x < grid_width - diag; ++grid_x) {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);

            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
            }
        }

        for (; grid_y < grid_height - diag; ++grid_y) {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);

            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
            }
        }

        for (; grid_x > diag; --grid_x) {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);

            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
            }
        }

        {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);
            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
            }
        }

        for (; grid_y > diag; --grid_y) {
            const auto cell_counts =
                count_cells_in(tissue, t_bbox, grid_x, grid_y, width, height);
            if (constraints_satisfied(cell_counts, species_constraints,
                                      mutant_constraints)) {
                return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
            }
        }
    }
    Rcpp::stop("No bounding box found!");
}

Logics::Variable TissueSimulation::get_var(const std::string &name) const
{

    if (name == "Time") {
        return Logics::Variable(sim_ptr->get_time_variable());
    }

    auto dot_pos = name.find('.');

    if (dot_pos == std::string::npos) {
        return Logics::Variable(sim_ptr->get_cardinality_variable(name));
    }

    const std::string event_name(name.substr(dot_pos + 1));
    const std::string species_name{name.substr(0, dot_pos)};

    const auto event_id = get_event_id(event_name, cell_event_names_inv);

    return sim_ptr->get_event_variable(species_name, event_id);
}
