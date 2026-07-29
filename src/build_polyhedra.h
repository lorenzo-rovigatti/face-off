#pragma once

#include <vector>
#include <tuple>
#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <iostream>
#include <cstdint>

// ============================================================================
// Triangle class: encapsulates the state of a single triangle in the capsid
// This class is designed for performance - it uses minimal storage and 
// provides direct access to state members for fast queries in the simulation.
// ============================================================================
class Triangle {
public:
    int num_inter_bonds;  // Number of bonds to other triangles (0-3)
    int num_polymer_bonds; // Number of bonds to polymer (0-1)
    bool is_remaining;         // Whether this triangle is still in the capsid
    
    // Constructor
    explicit Triangle(int initial_inter_bonds, int initial_polymer_bonds) 
        : num_inter_bonds(initial_inter_bonds), num_polymer_bonds(initial_polymer_bonds), is_remaining(true) {}
    
    // Get total number of bonds
    int total_bonds() const {
        return num_inter_bonds + num_polymer_bonds;
    }
};

using Face = std::tuple<int,int,int>;

using NeighborSideMap = std::vector<std::map<int, int>>;

std::vector<Face> icosahedron_faces() {
    return {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1},
    };
}

NeighborSideMap icosahedron_neighbor_sides() {
    std::vector<Face> faces = icosahedron_faces();
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edge_to_faces;

    for(uint32_t f_id = 0; f_id < faces.size(); ++f_id) {
        auto [v0, v1, v2] = faces[f_id];
        std::vector<int> vertices = {v0, v1, v2};

        for(int side = 0; side < 3; side++) {
            int a = vertices[side];
            int b = vertices[(side + 1) % 3];
            std::pair<int, int> edge = {std::min(a, b), std::max(a, b)};
            edge_to_faces[edge].push_back({(int)f_id, side});
        }
    }

    NeighborSideMap neighbor_sides(faces.size());
    for(const auto& [edge, attached] : edge_to_faces) {
        if(attached.size() == 2) {
            auto [i, side_i] = attached[0];
            auto [j, side_j] = attached[1];
            neighbor_sides[i][j] = side_i;
            neighbor_sides[j][i] = side_j;
        }
    }

    return neighbor_sides;
}

std::vector<std::vector<int>> icosahedron_neighbors() {
    std::vector<Face> faces = icosahedron_faces();
    std::map<std::pair<int, int>, std::vector<int>> edge_to_faces;
    
    for(uint32_t f_id = 0; f_id < faces.size(); ++f_id) {
        auto [v0, v1, v2] = faces[f_id];
        std::vector<int> vertices = {v0, v1, v2};
        
        for(int i = 0; i < 3; i++) {
            int a = vertices[i];
            int b = vertices[(i + 1) % 3];
            std::pair<int, int> edge = {std::min(a, b), std::max(a, b)};
            edge_to_faces[edge].push_back(f_id);
        }
    }

    std::map<int, std::set<int>> neighbors_set;
    for(uint32_t i = 0; i < faces.size(); i++) {
        neighbors_set[i] = std::set<int>();
    }

    for(const auto& [edge, attached] : edge_to_faces) {
        if(attached.size() == 2) {
            int i = attached[0];
            int j = attached[1];
            neighbors_set[i].insert(j);
            neighbors_set[j].insert(i);
        }
    }

    // Convert sets to sorted vectors
    std::vector<std::vector<int>> neighbors(faces.size());
    for(uint32_t i = 0; i < faces.size(); i++) {
        neighbors[i] = std::vector<int>(neighbors_set[i].begin(), neighbors_set[i].end());
        std::sort(neighbors[i].begin(), neighbors[i].end());
    }

    return neighbors;
}

std::vector<std::vector<int>> AAV_neighbors() {
    auto add_bond = [](std::vector<std::vector<int>>& neigh, int i, int j) {
        neigh[i].push_back(j);
        neigh[j].push_back(i);
    };

    const int n_pentamers = 12;
    const int units_per_pentamer = 5;
    const int N = n_pentamers * units_per_pentamer; // 60
    auto faces20 = icosahedron_faces();

    std::vector<std::vector<int>> neigh(N);

    // Build icosahedron vertex graph: 12 vertices, each degree 5
    std::vector<std::set<int>> ico_neigh(n_pentamers);

    for(const auto& face : faces20) {
        auto [a,b,c] = face;

        ico_neigh[a].insert(b);
        ico_neigh[b].insert(a);

        ico_neigh[b].insert(c);
        ico_neigh[c].insert(b);

        ico_neigh[c].insert(a);
        ico_neigh[a].insert(c);
    }

    // Convert icosahedron neighbour sets to sorted vectors
    std::vector<std::vector<int>> ico_order(n_pentamers);

    for(int v = 0; v < n_pentamers; v++) {
        ico_order[v] = std::vector<int>(
            ico_neigh[v].begin(),
            ico_neigh[v].end()
        );

        if(ico_order[v].size() != 5) {
            std::cerr << "Error: icosahedron vertex "
                      << v << " has degree "
                      << ico_order[v].size() << "\n";
        }
    }

    // Map: for pentamer p, which slot points toward neighbouring pentamer q?
    std::map<std::pair<int,int>, int> slot_of_neighbor;

    for(int p = 0; p < n_pentamers; p++) {
        for(int slot = 0; slot < 5; slot++) {
            int q = ico_order[p][slot];
            slot_of_neighbor[{p,q}] = slot;
        }
    }

    // 1. Intra-pentamer bonds: each pentamer is a 5-cycle
    for(int p = 0; p < n_pentamers; p++) {
        for(int s = 0; s < 5; s++) {
            int i = 5*p + s;
            int j = 5*p + ((s + 1) % 5);

            add_bond(neigh, i, j);
        }
    }

    // 2. Inter-pentamer bonds: one bond for each icosahedron edge
    for(int p = 0; p < n_pentamers; p++) {
        for(int q : ico_order[p]) {
            if(p < q) { // Avoid double counting edges
                int sp = slot_of_neighbor[{p,q}];
                int sq = slot_of_neighbor[{q,p}];

                int i = 5*p + sp;
                int j = 5*q + sq;

                add_bond(neigh, i, j);
            }
        }
    }

    return neigh;
}

NeighborSideMap AAV_neighbor_sides() {
    const int n_pentamers = 12;
    const int units_per_pentamer = 5;
    const int N = n_pentamers * units_per_pentamer;
    auto faces20 = icosahedron_faces();

    std::vector<std::set<int>> ico_neigh(n_pentamers);
    for(const auto& face : faces20) {
        auto [a, b, c] = face;

        ico_neigh[a].insert(b);
        ico_neigh[b].insert(a);

        ico_neigh[b].insert(c);
        ico_neigh[c].insert(b);

        ico_neigh[c].insert(a);
        ico_neigh[a].insert(c);
    }

    std::vector<std::vector<int>> ico_order(n_pentamers);
    for(int v = 0; v < n_pentamers; v++) {
        ico_order[v] = std::vector<int>(ico_neigh[v].begin(), ico_neigh[v].end());
    }

    std::map<std::pair<int, int>, int> slot_of_neighbor;
    for(int p = 0; p < n_pentamers; p++) {
        for(int slot = 0; slot < 5; slot++) {
            slot_of_neighbor[{p, ico_order[p][slot]}] = slot;
        }
    }

    NeighborSideMap neighbor_sides(N);

    for(int p = 0; p < n_pentamers; p++) {
        for(int s = 0; s < units_per_pentamer; s++) {
            int unit = 5 * p + s;
            int left = 5 * p + ((s + units_per_pentamer - 1) % units_per_pentamer);
            int right = 5 * p + ((s + 1) % units_per_pentamer);
            int outward_pentamer = ico_order[p][s];
            int outward_slot = slot_of_neighbor[{outward_pentamer, p}];
            int outward = 5 * outward_pentamer + outward_slot;

            neighbor_sides[unit][left] = 0;
            neighbor_sides[unit][right] = 1;
            neighbor_sides[unit][outward] = 2;
        }
    }

    return neighbor_sides;
}

std::vector<std::vector<int>> deltoidal_hexecontahedron_neighbors() {
    return {
        {1, 2, 4, 6},  // face 0
        {0, 3, 5, 7},  // face 1
        {0, 3, 4, 8},  // face 2
        {1, 2, 5, 9},  // face 3
        {0, 2, 11, 13},  // face 4
        {1, 3, 12, 14},  // face 5
        {0, 10, 11, 16},  // face 6
        {1, 10, 12, 17},  // face 7
        {2, 13, 15, 18},  // face 8
        {3, 14, 15, 19},  // face 9
        {6, 7, 20, 21},  // face 10
        {4, 6, 16, 22},  // face 11
        {5, 7, 17, 23},  // face 12
        {4, 8, 18, 24},  // face 13
        {5, 9, 19, 25},  // face 14
        {8, 9, 26, 27},  // face 15
        {6, 11, 20, 28},  // face 16
        {7, 12, 21, 29},  // face 17
        {8, 13, 26, 30},  // face 18
        {9, 14, 27, 31},  // face 19
        {10, 16, 21, 32},  // face 20
        {10, 17, 20, 33},  // face 21
        {11, 24, 28, 34},  // face 22
        {12, 25, 29, 35},  // face 23
        {13, 22, 30, 36},  // face 24
        {14, 23, 31, 37},  // face 25
        {15, 18, 27, 38},  // face 26
        {15, 19, 26, 39},  // face 27
        {16, 22, 34, 40},  // face 28
        {17, 23, 35, 41},  // face 29
        {18, 24, 36, 42},  // face 30
        {19, 25, 37, 43},  // face 31
        {20, 33, 40, 44},  // face 32
        {21, 32, 41, 44},  // face 33
        {22, 28, 36, 45},  // face 34
        {23, 29, 37, 46},  // face 35
        {24, 30, 34, 47},  // face 36
        {25, 31, 35, 48},  // face 37
        {26, 39, 42, 49},  // face 38
        {27, 38, 43, 49},  // face 39
        {28, 32, 45, 50},  // face 40
        {29, 33, 46, 51},  // face 41
        {30, 38, 47, 52},  // face 42
        {31, 39, 48, 53},  // face 43
        {32, 33, 50, 51},  // face 44
        {34, 40, 50, 54},  // face 45
        {35, 41, 51, 55},  // face 46
        {36, 42, 52, 54},  // face 47
        {37, 43, 53, 55},  // face 48
        {38, 39, 52, 53},  // face 49
        {40, 44, 45, 56},  // face 50
        {41, 44, 46, 57},  // face 51
        {42, 47, 49, 58},  // face 52
        {43, 48, 49, 59},  // face 53
        {45, 47, 56, 58},  // face 54
        {46, 48, 57, 59},  // face 55
        {50, 54, 57, 58},  // face 56
        {51, 55, 56, 59},  // face 57
        {52, 54, 56, 59},  // face 58
        {53, 55, 57, 58},  // face 59
    };
}

NeighborSideMap deltoidal_hexecontahedron_neighbor_sides() {
    auto neighbors = deltoidal_hexecontahedron_neighbors();
    NeighborSideMap neighbor_sides(neighbors.size());

    for(size_t i = 0; i < neighbors.size(); i++) {
        if(neighbors[i].size() != 4) {
            std::cerr << "Error: deltoidal face " << i
                      << " does not have 4 neighbors (found "
                      << neighbors[i].size() << ")\n";
        }

        for(size_t side = 0; side < neighbors[i].size(); side++) {
            int j = neighbors[i][side];
            neighbor_sides[i][j] = (int)side;
        }
    }

    return neighbor_sides;
}
