// Directed road graph in compressed sparse row form.
//
// Nodes are road junctions and dead ends; the points in between are kept
// only as geometry for drawing. Edge weights are travel times in tenths of
// a second (rounded up, so they never undercut the straight-line bound
// used by A*), and lengths are in metres.
#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

constexpr uint32_t kInf = std::numeric_limits<uint32_t>::max();
constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();

struct Graph {
    uint32_t n = 0;
    std::vector<int32_t> lat, lon;           // 1e-7 degrees
    std::vector<uint32_t> first_out;         // size n + 1
    std::vector<uint32_t> head;              // target of each edge
    std::vector<uint32_t> weight;            // deciseconds
    std::vector<float> length;               // metres
    std::vector<uint32_t> geom_first;        // size m + 1, into geom_lat/lon
    std::vector<int32_t> geom_lat, geom_lon; // intermediate points of each edge

    // reverse graph: for each node, its incoming edges
    std::vector<uint32_t> rfirst_in, rtail, redge;

    uint32_t m() const { return uint32_t(head.size()); }
    double max_speed_mps = 0;  // fastest edge, for the A* bound

    void build_reverse();
    void compute_max_speed();
    uint32_t find_edge(uint32_t u, uint32_t v, uint32_t w) const;  // kNone if absent

    void save(const std::string& path) const;
    static Graph load(const std::string& path);
};

// Simple edge list used while building.
struct EdgeIn {
    uint32_t u, v, w;
    float len;
    std::vector<int32_t> glat, glon;
};

Graph graph_from_edges(uint32_t n, std::vector<int32_t> lat, std::vector<int32_t> lon,
                       std::vector<EdgeIn>& edges);
