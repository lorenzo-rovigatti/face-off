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

// Structure to hold capsid parameters and precomputed rates
struct CapsidParameters {
    uint32_t N;
    std::vector<std::vector<int>> neighbors;
    double k_off;  // Per-bond breaking rate: nu * exp(-beta_eps)
    // Future extension: add k_on (per-bond re-forming rate) here

    CapsidParameters(const std::vector<std::vector<int>>& mneighbors, double beta_eps, double nu):
            N(mneighbors.size()),
            neighbors(mneighbors),
            k_off(nu * std::exp(-beta_eps))
    {}
};

// Structure to hold simulation results
struct SimulationResult {
    std::vector<double> time;
    double waiting_time;
    std::vector<int> n;            // Triangles remaining after each bond-break event
    std::vector<int> b;            // Bonds remaining after each bond-break event
    std::vector<int> event_bond_i; // Triangle i of the broken bond
    std::vector<int> event_bond_j; // Triangle j of the broken bond

    SimulationResult(uint32_t max_events) {
        time.reserve(max_events);
        n.reserve(max_events);
        b.reserve(max_events);
        event_bond_i.reserve(max_events);
        event_bond_j.reserve(max_events);
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

// Main Gillespie disassembly simulation: events are single bond breaks.
// Triangles that lose all bonds detach immediately as a consequence.
SimulationResult gillespie_disassembly(const CapsidParameters& params, unsigned int seed) {
    std::mt19937 gen(seed != 0 ? seed : std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Bond state: for each triangle, the set of neighbours it is currently bonded to.
    // To add re-forming bonds in the future, track broken bonds in a parallel set here.
    std::vector<std::set<int>> intact_bonds(params.N);
    std::vector<bool> is_remaining(params.N, true);

    for(uint32_t i = 0; i < params.N; i++) {
        for(int j : params.neighbors[i]) {
            intact_bonds[i].insert(j);
        }
    }

    int remaining_bonds = count_intact_bonds(intact_bonds);
    int remaining_triangles = params.N;

    // Each bond breaks at most once, so pre-allocate for that many events.
    SimulationResult result(remaining_bonds);

    double t = 0.0;

    // Gillespie loop: each event is the breaking of one bond.
    while(remaining_bonds > 0) {
        // All intact bonds break at the same rate k_off; total rate scales linearly.
        double K = remaining_bonds * params.k_off;

        // Sample waiting time to next bond-break event.
        double dt = -std::log(uniform(gen)) / K;
        t += dt;

        // Pick a bond uniformly at random (all rates are equal).
        int bond_idx = std::min((int)(uniform(gen) * remaining_bonds), remaining_bonds - 1);

        // Walk the intact-bond adjacency to find the bond_idx-th bond (i < j).
        int chosen_i = -1, chosen_j = -1;
        int count = 0;
        for(uint32_t i = 0; i < params.N && chosen_i < 0; i++) {
            if(!is_remaining[i]) continue;
            for(int j : intact_bonds[i]) {
                if(j > (int)i) {
                    if(count == bond_idx) {
                        chosen_i = (int)i;
                        chosen_j = j;
                        break;
                    }
                    count++;
                }
            }
        }

        // Break the chosen bond on both endpoints.
        intact_bonds[chosen_i].erase(chosen_j);
        intact_bonds[chosen_j].erase(chosen_i);
        remaining_bonds--;

        // Detach any triangle that now has zero intact bonds.
        for(int tri : {chosen_i, chosen_j}) {
            if(is_remaining[tri] && intact_bonds[tri].empty()) {
                is_remaining[tri] = false;
                remaining_triangles--;
            }
        }

        // Record the event.
        result.time.push_back(t);
        result.n.push_back(remaining_triangles);
        result.b.push_back(remaining_bonds);
        result.event_bond_i.push_back(chosen_i);
        result.event_bond_j.push_back(chosen_j);
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
    fmt::print("{:>12} {:>8} {:>8} {:>8} {:>8}\n", "Time", "N", "B", "BondI", "BondJ");

    for(size_t i = 0; i < result.time.size(); i++) {
        fmt::print("{:>12.6f} {:>8d} {:>8d} {:>8d} {:>8d}\n",
                    result.time[i],
                    result.n[i],
                    result.b[i],
                    result.event_bond_i[i],
                    result.event_bond_j[i]);
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
        ("temperature,T", "Temperature", cxxopts::value<double>())
        ("eps", "Epsilon parameter", cxxopts::value<double>())
        ("nu", "Attempt frequency", cxxopts::value<double>()->default_value("1.0"))
        ("trajectories", "Number of trajectories to average", cxxopts::value<unsigned int>()->default_value("1"))
        ("seed", "Random seed for reproducibility", cxxopts::value<unsigned int>()->default_value("42"))
        ("capsid", "Capsid type: 'aav' or 'icosahedron'", cxxopts::value<std::string>()->default_value("icosahedron"))
        ("output", "Output filename for averaged results", cxxopts::value<std::string>()->default_value("results.dat"))
        ("help", "Print help");
    
    double beta_eps;
    double nu;
    unsigned int n_trajectories;
    unsigned int seed;
    std::string capsid_type;
    std::vector<std::vector<int>> neighbors;
    std::string filename;
    
    try {
        auto result = options.parse(argc, argv);
        
        if(result.count("help")) {
            fmt::print("{}\n", options.help());
            return 0;
        }
        
        // Determine beta_eps
        if(result.count("beta_eps")) {
            beta_eps = result["beta_eps"].as<double>();
        }
        else if(result.count("T") && result.count("eps")) {
            beta_eps = result["eps"].as<double>() / result["T"].as<double>();
        }
        else {
            fmt::print("Error: Either --beta_eps or both --T and --eps must be specified\n");
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
    fmt::print("Parameters: beta_eps={}, nu={}, trajectories={}, seed={}, capsid={}\n", beta_eps, nu, n_trajectories, seed, capsid_type);
    fmt::print("\n");

    std::vector<SimulationResult> trajectories;
    CapsidParameters params(neighbors, beta_eps, nu);
    
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
