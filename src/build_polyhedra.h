#pragma once

#include <vector>
#include <tuple>
#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <iostream>

using Face = std::tuple<int,int,int>;

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