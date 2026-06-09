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
#include <fmt/core.h>
#include <fmt/format.h>

// Structure to hold some structural parameters of the capsid and precompute rates for different m values
struct CapsidParameters {
    int N;
    // Inter-triangle neighbor graph
    std::vector<std::vector<int>> neighbors;
    // Precomputed rates for different inter-bond counts
    std::vector<double> protein_rates;
    // Precomputed rates for different combined bond counts
    double polymer_rate;

    CapsidParameters(const std::vector<std::vector<int>> &mneighbors, double beta_eps_protein, double beta_eps_polymer, double nu): 
            N(mneighbors.size()), 
            neighbors(mneighbors) {

        // Precompute rates for inter-triangle bonds (0-3)
        int max_m = 0;
        for(int i = 0; i < N; i++) {
            max_m = std::max(max_m, static_cast<int>(neighbors[i].size()));
        }

        protein_rates.resize(max_m + 1);
        for(int m = 0; m <= max_m; m++) {
            protein_rates[m] = nu * std::exp(-beta_eps_protein * m);
        }
        
        polymer_rate = nu * std::exp(-beta_eps_polymer);
    }
};

// Structure to hold simulation results
struct SimulationResult {
    std::vector<double> time;
    double waiting_time;
    std::vector<int> n;                    // Number of remaining triangles
    std::vector<int> b_inter;              // Number of remaining inter-triangle bonds
    std::vector<int> b_polymer;            // Number of remaining polymer bonds
    std::vector<int> removed_triangle;     // Which triangle was removed
    std::vector<int> broken_inter_bonds;   // Number of inter-bonds broken in each event
    std::vector<int> broken_polymer_bonds; // Number of polymer bonds broken in each event

    SimulationResult(int N) {
        time.reserve(N);
        n.reserve(N);
        b_inter.reserve(N);
        b_polymer.reserve(N);
        removed_triangle.reserve(N);
        broken_inter_bonds.reserve(N);
        broken_polymer_bonds.reserve(N);
    }
};

// Structure to hold averaged results
struct AveragedResult {
    AveragedResult(int N_time_grid, double max_time) {
        // build the common time grid
        double dt = max_time / (N_time_grid - 1);
        time.resize(N_time_grid);
        for(int i = 0; i < N_time_grid; i++) {
            time[i] = dt * i;
        }

        // Average n and bonds at each step
        n_mean.resize(N_time_grid, 0.0);
        n_std.resize(N_time_grid, 0.0);
        b_inter_mean.resize(N_time_grid, 0.0);
        b_inter_std.resize(N_time_grid, 0.0);
        b_polymer_mean.resize(N_time_grid, 0.0);
        b_polymer_std.resize(N_time_grid, 0.0);
    }

    std::vector<double> time;
    double waiting_time_mean;
    double waiting_time_std;
    std::vector<double> n_mean;
    std::vector<double> n_std;
    std::vector<double> b_inter_mean;
    std::vector<double> b_inter_std;
    std::vector<double> b_polymer_mean;
    std::vector<double> b_polymer_std;
};

// Count intact inter-triangle bonds: slow function, use only for debugging or initial count
int count_intact_inter_bonds(const std::vector<Triangle>& triangles, 
                              const std::vector<std::vector<int>>& neighbors) {
    int b = 0;
    for(size_t i = 0; i < triangles.size(); i++) {
        if(triangles[i].num_inter_bonds > 0) {
            for(int j : neighbors[i]) {
                if(triangles[j].num_inter_bonds > 0 && j > (int)i) {
                    b++;
                }
            }
        }
    }
    return b;
}

// Count total intact polymer bonds: slow function, use only for debugging or initial count
int count_intact_polymer_bonds(const std::vector<Triangle>& triangles) {
    int b = 0;
    for(const auto& t : triangles) {
        b += t.num_polymer_bonds;
    }
    return b;
}

// Main Gillespie disassembly simulation
enum EventType : uint8_t {
    LOSE_INTER_BONDS = 0,   // Triangle detaches from capsid
    LOSE_POLYMER_BOND = 1   // Triangle detaches from polymer
};

// Compact event representation: which triangle, which type of event
struct CandidateEvent {
    int triangle_id;
    EventType event_type;
};

// Main Gillespie disassembly simulation
SimulationResult gillespie_disassembly(const CapsidParameters& params, double nu, unsigned int seed) {
    std::mt19937 gen(seed != 0 ? seed : std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Initialize triangles with 3 protein bonds and 1 polymer bond each
    std::vector<Triangle> triangles(params.N, Triangle(3, 1));

    double t = 0.0;
    SimulationResult result(params.N + 1);
    result.time.push_back(0.0);
    result.broken_inter_bonds.push_back(0);
    result.broken_polymer_bonds.push_back(0);
    
    int remaining_inter_bonds = count_intact_inter_bonds(triangles, params.neighbors);
    int remaining_polymer_bonds = count_intact_polymer_bonds(triangles);
    int fully_detached_triangles = 0;  // Count of triangles with no bonds at all

    // Storage for events (up to 2 events per triangle: lose inter-bonds or lose polymer-bond)
    std::vector<CandidateEvent> events;
    events.reserve(params.N * 2);
    std::vector<double> event_rates;
    event_rates.reserve(params.N * 2);

    // Loop continues while there are triangles with at least one bond (inter or polymer)
    while(fully_detached_triangles < params.N) {
        events.clear();
        event_rates.clear();

        // For each triangle, enumerate all possible competing events
        for(int i = 0; i < params.N; i++) {
            // Triangle must have at least one bond to have events
            bool has_inter_bonds = triangles[i].num_inter_bonds > 0;
            bool has_polymer_bonds = triangles[i].num_polymer_bonds > 0;
            
            // Event A: Triangle loses inter-bonds (detaches from capsid, but may stay via polymer)
            if(has_inter_bonds) {
                // Count intact inter-triangle bonds for this triangle
                // Only count neighbors that still have inter-bonds (in capsid)
                int m_inter = 0;
                for(int j : params.neighbors[i]) {
                    if(triangles[j].num_inter_bonds > 0) {
                        m_inter++;
                    }
                }
                
                double rate_lose_inter = params.protein_rates[m_inter];
                events.push_back({i, LOSE_INTER_BONDS});
                event_rates.push_back(rate_lose_inter);
            }
            
            // Event B: Triangle loses polymer-bond (stays in system via capsid or loses everything)
            if(has_polymer_bonds) {
                events.push_back({i, LOSE_POLYMER_BOND});
                event_rates.push_back(params.polymer_rate);
            }
        }

        // Calculate total rate (sum of all competing events)
        double K = std::accumulate(event_rates.begin(), event_rates.end(), 0.0);

        // Sample time to next event
        double u = uniform(gen);
        double dt = -std::log(u) / K;
        t += dt;

        // Sample which event fires using categorical distribution
        double u2 = uniform(gen);
        double cumsum = 0.0;
        size_t chosen_event_idx = 0;
        for(size_t i = 0; i < event_rates.size(); i++) {
            cumsum += event_rates[i] / K;  // normalize to get cumulative probabilities
            if(u2 <= cumsum) {
                chosen_event_idx = i;
                break;
            }
        }

        // Execute the chosen event
        const CandidateEvent& event = events[chosen_event_idx];
        int triangle_id = event.triangle_id;
        
        if(event.event_type == LOSE_INTER_BONDS) {
            // Triangle detaches from capsid structure
            // Count how many inter-bonds it had
            uint8_t m_inter_lost = 0;
            for(int j : params.neighbors[triangle_id]) {
                if(triangles[j].num_inter_bonds > 0) {
                    m_inter_lost++;
                }
            }
            
            remaining_inter_bonds -= m_inter_lost;
            triangles[triangle_id].num_inter_bonds = 0;
            
            result.broken_inter_bonds.push_back(m_inter_lost);
            result.broken_polymer_bonds.push_back(0);
        } 
        else if(event.event_type == LOSE_POLYMER_BOND) {
            // Triangle releases from polymer
            remaining_polymer_bonds -= triangles[triangle_id].num_polymer_bonds;
            triangles[triangle_id].num_polymer_bonds = 0;

            result.broken_inter_bonds.push_back(0);
            result.broken_polymer_bonds.push_back(1);
        }

        // Check if triangle is now fully detached (no bonds at all)
        if(triangles[triangle_id].total_bonds() == 0) {
            fully_detached_triangles++;
        }

        result.time.push_back(t);
        result.n.push_back(params.N - fully_detached_triangles);
        result.b_inter.push_back(remaining_inter_bonds);
        result.b_polymer.push_back(remaining_polymer_bonds);
        result.removed_triangle.push_back(triangle_id);
    }

    // Adjust times so that t=0 coincides with first event
    double t_wait = result.time[1];
    result.waiting_time = t_wait;
    for(size_t i = 1; i < result.time.size(); i++) {
        result.time[i] -= t_wait;
    }
    result.time.erase(result.time.begin());

    return result;
}

// Print simulation results
void print_results(const SimulationResult& result) {
    fmt::print("Waiting time: {:.6f}\n", result.waiting_time);
    fmt::print("\nSimulation timeline:\n");
    fmt::print("{:>12} {:>8} {:>9} {:>9} {:>12} {:>14} {:>14}\n", 
               "Time", "N", "B_inter", "B_poly", "Removed", "BrokenInter", "BrokenPoly");
    
    for(size_t i = 0; i < result.time.size(); i++) {
        fmt::print("{:>12.6f} {:>8d} {:>9d} {:>9d} {:>12d} {:>14d} {:>14d}\n",
                    result.time[i],
                    result.n[i],
                    result.b_inter[i],
                    result.b_polymer[i],
                    result.removed_triangle[i],
                    result.broken_inter_bonds[i],
                    result.broken_polymer_bonds[i]);
    }
}

// Compute averages over multiple trajectories
AveragedResult compute_averaged_results(const std::vector<SimulationResult>& trajectories, int N_time_grid) {
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
    std::vector<double> b_inter_sq_sum(N_time_grid, 0.0);
    std::vector<double> b_polymer_sq_sum(N_time_grid, 0.0);

    for(const auto& traj : trajectories) {
        int current_common_idx = 0;
        double common_time = averaged.time[current_common_idx];
         for(int i = 1; i < (int)traj.n.size(); i++) {
            double t = traj.time[i];
            int prev_idx = i - 1;
            while(common_time < t) {
                averaged.n_mean[current_common_idx] += traj.n[prev_idx];
                averaged.n_std[current_common_idx] += traj.n[prev_idx] * traj.n[prev_idx];
                averaged.b_inter_mean[current_common_idx] += traj.b_inter[prev_idx];
                averaged.b_inter_std[current_common_idx] += traj.b_inter[prev_idx] * traj.b_inter[prev_idx];
                averaged.b_polymer_mean[current_common_idx] += traj.b_polymer[prev_idx];
                averaged.b_polymer_std[current_common_idx] += traj.b_polymer[prev_idx] * traj.b_polymer[prev_idx];

                current_common_idx++;
                common_time = averaged.time[current_common_idx];
            } 
         }
    }

    // Compute mean and std
    for(int i = 0; i < N_time_grid; i++) {
        averaged.n_mean[i] /= n_traj;
        averaged.b_inter_mean[i] /= n_traj;
        averaged.b_polymer_mean[i] /= n_traj;
        
        double n_var = (averaged.n_std[i] / n_traj) - (averaged.n_mean[i] * averaged.n_mean[i]);
        double b_inter_var = (averaged.b_inter_std[i] / n_traj) - (averaged.b_inter_mean[i] * averaged.b_inter_mean[i]);
        double b_polymer_var = (averaged.b_polymer_std[i] / n_traj) - (averaged.b_polymer_mean[i] * averaged.b_polymer_mean[i]);
        
        averaged.n_std[i] = std::sqrt(std::max(0.0, n_var));
        averaged.b_inter_std[i] = std::sqrt(std::max(0.0, b_inter_var));
        averaged.b_polymer_std[i] = std::sqrt(std::max(0.0, b_polymer_var));
    }

    return averaged;
}

// Print averaged results
void print_averaged_results(const std::string& filename, const AveragedResult& averaged) {
    std::ofstream ofs(filename);
    
    ofs << fmt::format("# Waiting time: {:.6f} ± {:.6f}", averaged.waiting_time_mean, averaged.waiting_time_std) << std::endl;
    ofs << fmt::format("# {:>12} {:>16} {:>16} {:>16} {:>16} {:>18} {:>18}", 
                       "Time", "N_mean", "N_std", "B_inter_mean", "B_inter_std", "B_polymer_mean", "B_polymer_std") << std::endl;
    
    for(size_t i = 0; i < averaged.time.size(); i++) {
        ofs << fmt::format("{:>12.6f} {:>16.2f} {:>16.2f} {:>16.2f} {:>16.2f} {:>18.2f} {:>18.2f}", 
                averaged.time[i], averaged.n_mean[i], averaged.n_std[i], 
                averaged.b_inter_mean[i], averaged.b_inter_std[i],
                averaged.b_polymer_mean[i], averaged.b_polymer_std[i])
            << std::endl;
    }

    ofs.close();
}

int main(int argc, char* argv[]) {
    // Parse command line arguments with cxxopts
    cxxopts::Options options("face-off", "Gillespie disassembly simulation for viral capsids");
    
    options.add_options()
        ("temperature,T", "Temperature", cxxopts::value<double>())
        ("eps_protein", "Protein-protein interaction parameter", cxxopts::value<double>())
        ("eps_polymer", "Polymer-protein interaction parameter", cxxopts::value<double>())
        ("nu", "Attempt frequency", cxxopts::value<double>()->default_value("1.0"))
        ("trajectories", "Number of trajectories to average", cxxopts::value<unsigned int>()->default_value("1"))
        ("seed", "Random seed for reproducibility", cxxopts::value<unsigned int>()->default_value("42"))
        ("capsid", "Capsid type: 'aav' or 'icosahedron'", cxxopts::value<std::string>()->default_value("icosahedron"))
        ("output", "Output filename for averaged results", cxxopts::value<std::string>()->default_value("results.dat"))
        ("help", "Print help");
    
    double T;
    double beta_eps_protein;
    double beta_eps_polymer;
    double nu;
    unsigned int n_trajectories;
    unsigned int seed;
    std::string capsid_type;
    std::vector<std::vector<int>> neighbors;
    std::string filename = "results.dat";
    
    try {
        auto result = options.parse(argc, argv);
        
        if(result.count("help")) {
            fmt::print("{}\n", options.help());
            return 0;
        }
        
        // interaction parameters
        if(result.count("T") && result.count("eps_protein") && result.count("eps_polymer")) {
            T = result["T"].as<double>();
            beta_eps_protein = result["eps_protein"].as<double>() / T;
            beta_eps_polymer = result["eps_polymer"].as<double>() / T;
        }
        else {
            fmt::print("Error: --T, --eps_protein, and --eps_polymer must all be specified\n");
            fmt::print("{}\n", options.help());
            return 1;
        }
        
        nu = result["nu"].as<double>();
        n_trajectories = result["trajectories"].as<unsigned int>();
        seed = result["seed"].as<unsigned int>();
        capsid_type = result["capsid"].as<std::string>();
        filename = result["output"].as<std::string>();
        
        // Select capsid type
        if(capsid_type == "aav") {
            neighbors = AAV_neighbors();
        }
        else if(capsid_type == "icosahedron") {
            neighbors = icosahedron_neighbors();
        }
        else {
            fmt::print("Error: Unknown capsid type '{}'. Use 'aav' or 'icosahedron'\n", capsid_type);
            return 1;
        }
    }
    catch(const cxxopts::exceptions::exception& e) {
        fmt::print("Error parsing arguments: {}\n", e.what());
        fmt::print("{}\n", options.help());
        return 1;
    }

    fmt::print("Running Gillespie disassembly simulation...\n");
    fmt::print("Parameters: T={}, beta_eps_protein={}, beta_eps_polymer={}, nu={}, trajectories={}, seed={}, capsid={}\n", T, beta_eps_protein, beta_eps_polymer, nu, n_trajectories, seed, capsid_type);
    fmt::print("\n");

    std::vector<SimulationResult> trajectories;
    CapsidParameters params(neighbors, beta_eps_protein, beta_eps_polymer, nu);
    
    // Run multiple trajectories
    for(unsigned int i = 0; i < n_trajectories; i++) {
        unsigned int traj_seed = seed + i;
        auto result = gillespie_disassembly(params, nu, traj_seed);
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
