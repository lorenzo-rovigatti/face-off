#include "build_polyhedra.h"

#include <cxxopts.hpp>

#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <fmt/core.h>
#include <fmt/format.h>

// Structure to hold capsid parameters and precomputed rates
struct CapsidParameters {
    uint32_t N;
    std::vector<std::vector<int>> neighbors;
    NeighborSideMap neighbor_sides;
    std::map<std::pair<int, int>, double> k_off_by_bond;
    double k_on;   // Per-bond re-forming rate (T-independent, same units as nu; 0 = irreversible)

    CapsidParameters(const std::vector<std::vector<int>>& mneighbors,
                     const NeighborSideMap& mneighbor_sides,
                     const std::map<std::pair<int, int>, double>& beta_eps_by_pair,
                     double nu,
                     double k_on_param = 0.0):
            N(mneighbors.size()),
            neighbors(mneighbors),
            neighbor_sides(mneighbor_sides),
            k_on(k_on_param)
    {
        if(neighbor_sides.size() != neighbors.size()) {
            throw std::invalid_argument("Neighbor-side map size does not match topology size");
        }

        for(uint32_t i = 0; i < N; i++) {
            if(neighbor_sides[i].size() != neighbors[i].size()) {
                throw std::invalid_argument("Neighbor-side map does not cover all bonds");
            }

            for(int j : neighbors[i]) {
                if(j <= (int)i) continue;

                auto side_i_it = neighbor_sides[i].find(j);
                auto side_j_it = neighbor_sides[j].find(i);
                if(side_i_it == neighbor_sides[i].end() || side_j_it == neighbor_sides[j].end()) {
                    throw std::invalid_argument("Missing side assignment for bond");
                }

                int side_a = std::min(side_i_it->second, side_j_it->second);
                int side_b = std::max(side_i_it->second, side_j_it->second);
                auto pair_it = beta_eps_by_pair.find({side_a, side_b});
                if(pair_it == beta_eps_by_pair.end()) {
                    throw std::invalid_argument("Missing beta_eps value for bond type " + std::to_string(side_a) + "-" + std::to_string(side_b));
                }

                double beta_eps = pair_it->second;
                k_off_by_bond[{(int)i, j}] = nu * std::exp(-beta_eps);
            }
        }
    }

    double bond_break_rate(int i, int j) const {
        if(i > j) std::swap(i, j);
        return k_off_by_bond.at({i, j});
    }
};

std::string trim_copy(const std::string& text) {
    size_t first = text.find_first_not_of(" \t");
    if(first == std::string::npos) return "";
    size_t last = text.find_last_not_of(" \t");
    return text.substr(first, last - first + 1);
}

std::map<std::pair<int, int>, double> parse_pair_values(const std::string& text, const std::string& option_name) {
    std::map<std::pair<int, int>, double> values;
    std::stringstream stream(text);
    std::string token;

    while(std::getline(stream, token, ',')) {
        std::string entry = trim_copy(token);
        if(entry.empty()) {
            throw std::runtime_error("Option --" + option_name + " contains an empty assignment");
        }

        size_t colon_pos = entry.find(':');
        if(colon_pos == std::string::npos) {
            throw std::runtime_error("Option --" + option_name + " must use entries formatted as i-j:value");
        }

        std::string pair_text = trim_copy(entry.substr(0, colon_pos));
        std::string value_text = trim_copy(entry.substr(colon_pos + 1));

        size_t dash_pos = pair_text.find('-');
        if(dash_pos == std::string::npos) {
            throw std::runtime_error("Option --" + option_name + " must use side pairs formatted as i-j");
        }

        int side_i = std::stoi(trim_copy(pair_text.substr(0, dash_pos)));
        int side_j = std::stoi(trim_copy(pair_text.substr(dash_pos + 1)));

        if(side_i < 0 || side_i > 2 || side_j < 0 || side_j > 2) {
            throw std::runtime_error("Option --" + option_name + " side indices must be in [0,2]");
        }

        int side_a = std::min(side_i, side_j);
        int side_b = std::max(side_i, side_j);
        values[{side_a, side_b}] = std::stod(value_text);
    }

    if(values.empty()) {
        throw std::runtime_error("Option --" + option_name + " must contain at least one i-j:value assignment");
    }

    return values;
}

std::set<std::pair<int, int>> required_bond_types(const std::vector<std::vector<int>>& neighbors,
                                                  const NeighborSideMap& neighbor_sides) {
    std::set<std::pair<int, int>> required;
    for(size_t i = 0; i < neighbors.size(); i++) {
        for(int j : neighbors[i]) {
            if(j <= (int)i) continue;

            int side_i = neighbor_sides[i].at(j);
            int side_j = neighbor_sides[j].at((int)i);
            int side_a = std::min(side_i, side_j);
            int side_b = std::max(side_i, side_j);
            required.insert({side_a, side_b});
        }
    }
    return required;
}

std::string join_missing_bond_types(const std::set<std::pair<int, int>>& missing) {
    std::string out;
    bool first = true;
    for(const auto& [a, b] : missing) {
        if(!first) out += ", ";
        out += fmt::format("{}-{}", a, b);
        first = false;
    }
    return out;
}

// Structure to hold simulation results
struct SimulationResult {
    std::vector<double> time;
    double waiting_time;
    std::vector<int> n;                      // Triangles remaining after each event
    std::vector<int> b;                      // Bonds remaining after each event
    std::vector<int> event_bond_i;           // Triangle i of the affected bond
    std::vector<int> event_bond_j;           // Triangle j of the affected bond
    std::vector<bool> event_is_formation;    // true = bond formed, false = bond broken

    SimulationResult(uint32_t max_events) {
        time.reserve(max_events);
        n.reserve(max_events);
        b.reserve(max_events);
        event_bond_i.reserve(max_events);
        event_bond_j.reserve(max_events);
        event_is_formation.reserve(max_events);
    }
};

// Structure to hold averaged results
struct AveragedResult {
        AveragedResult(uint32_t N_time_grid, double max_time) {
            // build the common time grid
            double dt = max_time / (N_time_grid - 1);
            time.resize(N_time_grid);
            for(uint32_t i = 0; i < N_time_grid; i++) {
                time[i] = dt * i;
            }

            // Average n and bonds at each step
            n_mean.resize(N_time_grid, 0.0);
            n_std.resize(N_time_grid, 0.0);
            b_mean.resize(N_time_grid, 0.0);
            b_std.resize(N_time_grid, 0.0);
        }

    std::vector<double> time;
    double waiting_time_mean;
    double waiting_time_std;
    std::vector<double> n_mean;
    std::vector<double> n_std;
    std::vector<double> b_mean;
    std::vector<double> b_std;
};

// Count intact bonds from the bond-state sets (use i < j to avoid double-counting)
int count_intact_bonds(const std::vector<std::set<int>>& intact_bonds) {
    int b = 0;
    for(size_t i = 0; i < intact_bonds.size(); i++) {
        for(int j : intact_bonds[i]) {
            if(j > (int)i) b++;
        }
    }
    return b;
}

// Main Gillespie disassembly simulation: events are single bond breaks or re-formations.
// A broken bond can re-form (at rate k_on) only if both endpoint triangles still have at
// least one intact bond (i.e., neither has fully detached from the capsid).
// Triangles that lose all bonds detach permanently.
SimulationResult gillespie_disassembly(const CapsidParameters& params, unsigned int seed) {
    std::mt19937 gen(seed != 0 ? seed : std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<std::set<int>> intact_bonds(params.N);
    std::vector<bool> is_remaining(params.N, true);

    for(uint32_t i = 0; i < params.N; i++) {
        for(int j : params.neighbors[i]) {
            intact_bonds[i].insert(j);
        }
    }

    int remaining_bonds = count_intact_bonds(intact_bonds);
    int remaining_triangles = params.N;

    // With reversible bonds the event count is unbounded; reserve generously.
    SimulationResult result(remaining_bonds * 4);

    // Indexed set of re-formable bonds: O(1) random access, O(log n) insert/remove.
    // Invariant: contains bond (i,j) with i<j iff it is an original neighbour pair,
    // currently broken, and both endpoints are still in the capsid.
    //
    // Remove uses swap-with-last so the vector stays packed.
    std::vector<std::pair<int,int>> formable_vec;
    std::map<std::pair<int,int>, int> formable_idx;

    const auto formable_insert = [&](int i, int j) {
        if(i > j) std::swap(i, j);
        auto key = std::make_pair(i, j);
        formable_idx[key] = (int)formable_vec.size();
        formable_vec.push_back(key);
    };

    const auto formable_remove = [&](int i, int j) {
        if(i > j) std::swap(i, j);
        auto key = std::make_pair(i, j);
        auto it = formable_idx.find(key);
        if(it == formable_idx.end()) return; // not in set (bond may still be intact)
        int idx = it->second;
        int last = (int)formable_vec.size() - 1;
        if(idx != last) {
            // Swap with the last element and fix its index entry.
            formable_vec[idx] = formable_vec[last];
            formable_idx[formable_vec[idx]] = idx;
        }
        formable_vec.pop_back();
        formable_idx.erase(key);
    };

    double t = 0.0;

    // Gillespie loop: each event is either a bond break (k_off) or a bond re-formation (k_on).
    while(remaining_bonds > 0) {
        int n_formable = (int)formable_vec.size();
        double K_break = 0.0;
        for(uint32_t i = 0; i < params.N; i++) {
            if(!is_remaining[i]) continue;
            for(int j : intact_bonds[i]) {
                if(j > (int)i) {
                    K_break += params.bond_break_rate((int)i, j);
                }
            }
        }
        double K_form  = n_formable       * params.k_on;
        double K       = K_break + K_form;

        // Sample waiting time.
        double dt = -std::log(uniform(gen)) / K;
        t += dt;

        // Decide event type, then pick uniformly within that class.
        bool is_formation = (uniform(gen) * K >= K_break);

        int chosen_i = -1, chosen_j = -1;

        if(!is_formation) {
            // --- Break event ---
            double threshold = uniform(gen) * K_break;
            double cumulative = 0.0;
            for(uint32_t i = 0; i < params.N && chosen_i < 0; i++) {
                if(!is_remaining[i]) continue;
                for(int j : intact_bonds[i]) {
                    if(j > (int)i) {
                        cumulative += params.bond_break_rate((int)i, j);
                        if(cumulative >= threshold) { chosen_i = (int)i; chosen_j = j; break; }
                    }
                }
            }

            if(chosen_i < 0) {
                throw std::runtime_error("Failed to select a bond-breaking event");
            }

            intact_bonds[chosen_i].erase(chosen_j);
            intact_bonds[chosen_j].erase(chosen_i);
            remaining_bonds--;

            // Check each endpoint for detachment and update the formable set accordingly.
            bool i_detached = false, j_detached = false;
            if(intact_bonds[chosen_i].empty()) {
                is_remaining[chosen_i] = false;
                remaining_triangles--;
                i_detached = true;
                // All formable bonds that had chosen_i as an endpoint are now invalid.
                for(int k : params.neighbors[chosen_i])
                    formable_remove(chosen_i, k);
            }
            if(intact_bonds[chosen_j].empty()) {
                is_remaining[chosen_j] = false;
                remaining_triangles--;
                j_detached = true;
                for(int k : params.neighbors[chosen_j])
                    formable_remove(chosen_j, k);
            }
            // The broken bond itself becomes re-formable only if both ends are still attached.
            if(!i_detached && !j_detached)
                formable_insert(chosen_i, chosen_j);

        } 
        else {
            // --- Formation event ---
            int bond_idx = std::min((int)(uniform(gen) * n_formable), n_formable - 1);
            auto [ci, cj] = formable_vec[bond_idx];
            chosen_i = ci; chosen_j = cj;
            formable_remove(chosen_i, chosen_j);
            intact_bonds[chosen_i].insert(chosen_j);
            intact_bonds[chosen_j].insert(chosen_i);
            remaining_bonds++;
        }

        // Record the event.
        result.time.push_back(t);
        result.n.push_back(remaining_triangles);
        result.b.push_back(remaining_bonds);
        result.event_bond_i.push_back(chosen_i);
        result.event_bond_j.push_back(chosen_j);
        result.event_is_formation.push_back(is_formation);
    }

    // Shift times so that t=0 coincides with the first triangle detachment.
    double t_wait = result.time.back(); // fallback: no triangle ever detached
    for(size_t i = 0; i < result.n.size(); i++) {
        if(result.n[i] < (int)params.N) {
            t_wait = result.time[i];
            break;
        }
    }
    result.waiting_time = t_wait;
    for(double& ti : result.time) {
        ti -= t_wait;
    }

    return result;
}

// Print simulation results
void print_results(const SimulationResult& result) {
    fmt::print("Waiting time: {:.6f}\n", result.waiting_time);
    fmt::print("\nSimulation timeline:\n");
    fmt::print("{:>12} {:>8} {:>8} {:>8} {:>8} {:>8}\n", "Time", "N", "B", "BondI", "BondJ", "Event");

    for(size_t i = 0; i < result.time.size(); i++) {
        fmt::print("{:>12.6f} {:>8d} {:>8d} {:>8d} {:>8d} {:>8}\n",
                    result.time[i],
                    result.n[i],
                    result.b[i],
                    result.event_bond_i[i],
                    result.event_bond_j[i],
                    result.event_is_formation[i] ? "form" : "break");
    }
}

// Compute averages over multiple trajectories
AveragedResult compute_averaged_results(const std::vector<SimulationResult>& trajectories, uint32_t N_time_grid) {
    if(trajectories.empty()) {
        return AveragedResult(N_time_grid, 0.0);
    }

    size_t n_traj = trajectories.size();
    
    // Find maximum time across all trajectories to define common time grid
    double max_time = 0.0;
    for(const auto& traj : trajectories) {
        max_time = std::max(max_time, traj.time.back());
    }

    AveragedResult averaged(N_time_grid, max_time);

    // build the common time grid
    double dt = max_time / (N_time_grid - 1);
    averaged.time.resize(N_time_grid);
    for(uint32_t i = 0; i < N_time_grid; i++) {
        averaged.time[i] = dt * i;
    }
    
    // Average waiting times
    double waiting_time_sum = 0.0;
    double waiting_time_sq_sum = 0.0;
    for(const auto& traj : trajectories) {
        waiting_time_sum += traj.waiting_time;
        waiting_time_sq_sum += traj.waiting_time * traj.waiting_time;
    }
    averaged.waiting_time_mean = waiting_time_sum / n_traj;
    double waiting_time_var = (waiting_time_sq_sum / n_traj) - (averaged.waiting_time_mean * averaged.waiting_time_mean);
    averaged.waiting_time_std = std::sqrt(std::max(0.0, waiting_time_var));

    std::vector<double> n_sq_sum(N_time_grid, 0.0);
    std::vector<double> b_sq_sum(N_time_grid, 0.0);

    for(const auto& traj : trajectories) {
        uint32_t current_common_idx = 0;
        double common_time = averaged.time[current_common_idx];
         for(uint32_t i = 1; i < traj.n.size(); i++) {
            double t = traj.time[i];
            uint32_t prev_idx = i - 1;
            while(common_time < t) {
                averaged.n_mean[current_common_idx] += traj.n[prev_idx];
                averaged.n_std[current_common_idx] += traj.n[prev_idx] * traj.n[prev_idx];
                averaged.b_mean[current_common_idx] += traj.b[prev_idx];
                averaged.b_std[current_common_idx] += traj.b[prev_idx] * traj.b[prev_idx];

                current_common_idx++;
                common_time = averaged.time[current_common_idx];
            } 
         }
    }

    // Compute mean and std
    for(size_t i = 0; i < N_time_grid; i++) {
        averaged.n_mean[i] /= n_traj;
        averaged.b_mean[i] /= n_traj;
        
        double n_var = (averaged.n_std[i] / n_traj) - (averaged.n_mean[i] * averaged.n_mean[i]);
        double b_var = (averaged.b_std[i] / n_traj) - (averaged.b_mean[i] * averaged.b_mean[i]);
        
        averaged.n_std[i] = std::sqrt(std::max(0.0, n_var));
        averaged.b_std[i] = std::sqrt(std::max(0.0, b_var));
    }

    return averaged;
}

// Print averaged results
void print_averaged_results(const std::string& filename, const AveragedResult& averaged) {
    std::ofstream ofs(filename);
    
    ofs << fmt::format("# Waiting time: {:.6f} ± {:.6f}", averaged.waiting_time_mean, averaged.waiting_time_std) << std::endl;
    ofs << fmt::format("# {:>12} {:>16} {:>16} {:>16} {:>16}", "Time", "N_mean", "N_std", "B_mean", "B_std") << std::endl;
    
    for(size_t i = 0; i < averaged.time.size(); i++) {
        ofs << fmt::format("{:>12.6f} {:>16.2f} {:>16.2f} {:>16.2f} {:>16.2f}", 
                averaged.time[i], averaged.n_mean[i], averaged.n_std[i], averaged.b_mean[i], averaged.b_std[i])
            << std::endl;
    }

    ofs.close();
}

int main(int argc, char* argv[]) {
    // Parse command line arguments with cxxopts
    cxxopts::Options options("face-off", "Gillespie disassembly simulation for viral capsids");
    
    options.add_options()
        ("beta_eps", "Energy parameter (alternative to specifying beta and eps separately)", cxxopts::value<double>())
        ("beta_eps_bonds", "Bond-type beta*epsilon assignments as side pairs (i-j:value), e.g. 0-1:1.0,2-2:1.3", cxxopts::value<std::string>())
        ("temperature,T", "Temperature", cxxopts::value<double>())
        ("eps", "Epsilon parameter", cxxopts::value<double>())
        ("eps_bonds", "Bond-type epsilon assignments as side pairs (i-j:value); requires --T", cxxopts::value<std::string>())
        ("nu", "Attempt frequency", cxxopts::value<double>()->default_value("1.0"))
        ("k_on", "Bond re-formation rate constant (same units as nu; T-independent; 0 = irreversible)", cxxopts::value<double>()->default_value("0.0"))
        ("trajectories", "Number of trajectories to average", cxxopts::value<unsigned int>()->default_value("1"))
        ("seed", "Random seed for reproducibility", cxxopts::value<unsigned int>()->default_value("42"))
        ("output", "Output filename for averaged results", cxxopts::value<std::string>()->default_value("deltoid_results.dat"))
        ("help", "Print help");
    
    double beta_eps_global = 0.0;
    std::map<std::pair<int, int>, double> beta_eps_by_pair;
    double nu;
    double k_on;
    unsigned int n_trajectories;
    unsigned int seed;
    std::vector<std::vector<int>> neighbors;
    NeighborSideMap neighbor_sides;
    std::string filename;
    
    try {
        auto result = options.parse(argc, argv);
        
        if(result.count("help")) {
            fmt::print("{}\n", options.help());
            return 0;
        }
        
        bool has_global_beta_eps = result.count("beta_eps");
        bool has_global_eps = result.count("T") && result.count("eps");
        bool has_pair_beta_eps = result.count("beta_eps_bonds");
        bool has_pair_eps = result.count("eps_bonds");

        if((has_global_beta_eps || has_global_eps) && (has_pair_beta_eps || has_pair_eps)) {
            fmt::print("Error: Specify either global bond strength options or per-side bond strength options, not both\n");
            return 1;
        }

        if(has_pair_beta_eps) {
            beta_eps_by_pair = parse_pair_values(result["beta_eps_bonds"].as<std::string>(), "beta_eps_bonds");
        }
        else if(has_pair_eps) {
            if(!result.count("T")) {
                fmt::print("Error: --eps_bonds requires --T\n");
                return 1;
            }
            auto eps_by_pair = parse_pair_values(result["eps_bonds"].as<std::string>(), "eps_bonds");
            double temperature = result["T"].as<double>();
            for(const auto& [pair_key, eps_value] : eps_by_pair) {
                beta_eps_by_pair[pair_key] = eps_value / temperature;
            }
        }
        else if(has_global_beta_eps) {
            beta_eps_global = result["beta_eps"].as<double>();
        }
        else if(has_global_eps) {
            beta_eps_global = result["eps"].as<double>() / result["T"].as<double>();
        }
        else {
            fmt::print("Error: Specify either --beta_eps_bonds, --eps_bonds with --T, --beta_eps, or both --T and --eps\n");
            return 1;
        }
        
        nu = result["nu"].as<double>();
        k_on = result["k_on"].as<double>();
        n_trajectories = result["trajectories"].as<unsigned int>();
        seed = result["seed"].as<unsigned int>();
        filename = result["output"].as<std::string>();
        
        neighbors = deltoidal_hexecontahedron_neighbors();
        neighbor_sides = deltoidal_hexecontahedron_neighbor_sides();

        auto required_pairs = required_bond_types(neighbors, neighbor_sides);
        if(beta_eps_by_pair.empty()) {
            for(const auto& pair_key : required_pairs) {
                beta_eps_by_pair[pair_key] = beta_eps_global;
            }
        }
        else {
            std::set<std::pair<int, int>> missing;
            for(const auto& pair_key : required_pairs) {
                if(!beta_eps_by_pair.count(pair_key)) {
                    missing.insert(pair_key);
                }
            }
            if(!missing.empty()) {
                fmt::print("Error: Missing bond-type values for: {}\n", join_missing_bond_types(missing));
                return 1;
            }
        }
    }
    catch(const cxxopts::exceptions::exception& e) {
        fmt::print("Error parsing arguments: {}\n", e.what());
        fmt::print("{}\n", options.help());
        return 1;
    }

    fmt::print("Running Gillespie disassembly simulation...\n");
    std::string beta_eps_desc;
    bool first_pair = true;
    for(const auto& [pair_key, value] : beta_eps_by_pair) {
        if(!first_pair) beta_eps_desc += ",";
        beta_eps_desc += fmt::format("{}-{}:{}", pair_key.first, pair_key.second, value);
        first_pair = false;
    }
    fmt::print("Parameters: beta_eps_bonds={}, nu={}, k_on={}, trajectories={}, seed={}\n", beta_eps_desc, nu, k_on, n_trajectories, seed);
    fmt::print("\n");

    std::vector<SimulationResult> trajectories;
    CapsidParameters params(neighbors, neighbor_sides, beta_eps_by_pair, nu, k_on);
    
    // Run multiple trajectories
    for(unsigned int i = 0; i < n_trajectories; i++) {
        unsigned int traj_seed = seed + i;
        auto result = gillespie_disassembly(params, traj_seed);
        trajectories.push_back(result);
        
        if(n_trajectories <= 3) {
            fmt::print("Trajectory {}:\n", (i + 1));
            print_results(result);
            fmt::print("\n");
        } 
        else if((i + 1) % std::max(1u, n_trajectories / 10) == 0) {
            fmt::print("Completed {}/{} trajectories...\n", (i + 1), n_trajectories);
        }
    }

    if(n_trajectories > 1) {
        auto averaged = compute_averaged_results(trajectories, 1000);
        print_averaged_results(filename, averaged);
        fmt::print("\n=== AVERAGES PRINTED TO {} ===\n", filename);
    }

    return 0;
}
