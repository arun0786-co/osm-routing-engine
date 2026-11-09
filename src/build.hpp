#pragma once

#include <string>

#include "graph.hpp"
#include "pbf.hpp"

struct RoadClass {
    bool routable = false;
    double speed_kmh = 0;
    bool forward = true, backward = true;
};

// Car profile: which ways are drivable, how fast, and in which direction.
RoadClass classify(const osm::Tags& tags);

struct BuildStats {
    size_t ways = 0, refs = 0, osm_nodes = 0, missing_coords = 0;
    uint32_t junctions = 0, edges_before = 0;
    double seconds_ways = 0, seconds_nodes = 0, seconds_graph = 0;
};

Graph build_graph(const std::string& pbf_path, int threads, BuildStats* stats = nullptr);
