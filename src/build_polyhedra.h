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

std::vector<std::vector<int>> deltoidal_hexecontahedron_neighbors() {
    return {
        {2, 4, 6, 1}, // face 0
        {5, 3, 0, 7}, // face 1
        {4, 0, 3, 8}, // face 2
        {1, 5, 9, 2}, // face 3
        {0, 2, 13, 11}, // face 4
        {3, 1, 12, 14}, // face 5
        {11, 16, 10, 0}, // face 6
        {17, 12, 1, 10}, // face 7
        {18, 13, 2, 15}, // face 8
        {14, 19, 15, 3}, // face 9
        {20, 21, 7, 6}, // face 10
        {16, 6, 4, 22}, // face 11
        {7, 17, 23, 5}, // face 12
        {8, 18, 24, 4}, // face 13
        {19, 9, 5, 25}, // face 14
        {27, 26, 8, 9}, // face 15
        {6, 11, 28, 20}, // face 16
        {12, 7, 21, 29}, // face 17
        {13, 8, 26, 30}, // face 18
        {9, 14, 31, 27}, // face 19
        {21, 10, 16, 32}, // face 20
        {10, 20, 33, 17}, // face 21
        {34, 28, 11, 24}, // face 22
        {29, 35, 25, 12}, // face 23
        {30, 36, 22, 13}, // face 24
        {37, 31, 14, 23}, // face 25
        {15, 27, 38, 18}, // face 26
        {26, 15, 19, 39}, // face 27
        {22, 34, 40, 16}, // face 28
        {35, 23, 17, 41}, // face 29
        {36, 24, 18, 42}, // face 30
        {25, 37, 43, 19}, // face 31
        {44, 33, 20, 40}, // face 32
        {32, 44, 41, 21}, // face 33
        {28, 22, 36, 45}, // face 34
        {23, 29, 46, 37}, // face 35
        {24, 30, 47, 34}, // face 36
        {31, 25, 35, 48}, // face 37
        {39, 49, 42, 26}, // face 38
        {49, 38, 27, 43}, // face 39
        {45, 50, 32, 28}, // face 40
        {51, 46, 29, 33}, // face 41
        {52, 47, 30, 38}, // face 42
        {48, 53, 39, 31}, // face 43
        {33, 32, 50, 51}, // face 44
        {50, 40, 34, 54}, // face 45
        {41, 51, 55, 35}, // face 46
        {42, 52, 54, 36}, // face 47
        {53, 43, 37, 55}, // face 48
        {38, 39, 53, 52}, // face 49
        {40, 45, 56, 44}, // face 50
        {46, 41, 44, 57}, // face 51
        {47, 42, 49, 58}, // face 52
        {43, 48, 59, 49}, // face 53
        {58, 56, 45, 47}, // face 54
        {57, 59, 48, 46}, // face 55
        {54, 58, 57, 50}, // face 56
        {59, 55, 51, 56}, // face 57
        {56, 54, 52, 59}, // face 58
        {55, 57, 58, 53}, // face 59
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
