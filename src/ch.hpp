// Contraction hierarchies (Geisberger, Sanders, Schultes, Delling 2008).
//
// Preprocessing removes nodes one at a time in order of importance. When
// node v is removed, every path u -> v -> x that might be a shortest path
// is replaced by a shortcut u -> x, unless a local search finds a path from
// u to x avoiding v that is at least as short (a witness). Every shortest
// path then exists as a path that first only climbs to more important
// nodes and then only descends, so a query is two small upward searches.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "graph.hpp"
#include "search.hpp"

struct CHEdge {
    uint32_t to;
    uint32_t w;
    uint32_t mid;  // contracted node the shortcut skips, kNone for a road edge
};

struct CHParams {
    uint32_t witness_settle_sim = 60;    // settle limit while estimating priorities
    uint32_t witness_settle = 600;       // settle limit while contracting
    bool verbose = true;
};

struct CHStats {
    uint64_t shortcuts = 0;
    double seconds = 0;
    uint32_t max_level = 0;
};

struct CH {
    uint32_t n = 0;
    std::vector<uint32_t> rank;          // position in the contraction order
    // upward edges u -> v (rank v > rank u), searched forward from the source
    std::vector<uint32_t> up_first;
    std::vector<CHEdge> up;
    // edges u -> v stored at v when rank u > rank v; the backward search
    // from the target walks them from v to u
    std::vector<uint32_t> dn_first;
    std::vector<CHEdge> dn;

    static CH build(const Graph& g, const CHParams& p, CHStats* st = nullptr);
    void save(const std::string& path) const;
    static CH load(const std::string& path);
};

class CHQuery {
public:
    CHQuery(const CH& ch, const Graph& g);

    uint32_t run(uint32_t s, uint32_t t, QueryStats* st = nullptr);
    // road edges of the last query's shortest path, in order
    std::vector<uint32_t> edge_path() const;

    // distance table, sources x targets, row-major (bucket algorithm,
    // Knopp et al. 2007)
    std::vector<uint32_t> table(const std::vector<uint32_t>& sources,
                                const std::vector<uint32_t>& targets);

private:
    void upward(uint32_t s, Labels& lab, LazyHeap& heap, bool forward,
                std::vector<uint32_t>* settled);
    bool stalled(uint32_t v, uint32_t d, const Labels& lab, bool forward) const;
    void unpack(uint32_t u, const CHEdge& e, std::vector<uint32_t>& out) const;

    const CH& ch_;
    const Graph& g_;
    Labels f_, b_;
    LazyHeap hf_, hb_;
    uint32_t meet_ = kNone, s_ = kNone, t_ = kNone;
};
