// Baseline shortest-path searches on the plain road graph.
//
// All searches reuse their arrays between queries: a per-node timestamp
// says whether the stored distance belongs to the current query, so a
// query costs time proportional to what it touches, not to the graph.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <random>
#include <vector>

#include "geo.hpp"
#include "graph.hpp"

struct QueryStats {
    uint64_t settled = 0;
    uint64_t relaxed = 0;
};

// Min-heap of (key, node) packed into one 64-bit word. Stale entries are
// skipped when popped instead of being decreased in place.
class LazyHeap {
public:
    void clear() { h_.clear(); }
    bool empty() const { return h_.empty(); }
    void push(uint32_t key, uint32_t v) {
        h_.push_back((uint64_t(key) << 32) | v);
        std::push_heap(h_.begin(), h_.end(), std::greater<>());
    }
    uint32_t top_key() const { return uint32_t(h_.front() >> 32); }
    std::pair<uint32_t, uint32_t> pop() {
        std::pop_heap(h_.begin(), h_.end(), std::greater<>());
        uint64_t x = h_.back();
        h_.pop_back();
        return {uint32_t(x >> 32), uint32_t(x)};
    }

private:
    std::vector<uint64_t> h_;
};

// Distances with O(1) reset.
class Labels {
public:
    explicit Labels(uint32_t n) : dist_(n, kInf), parent_(n, kNone), stamp_(n, 0) {}
    void reset() {
        if (++now_ == 0) {  // wrapped: clear for real
            std::fill(stamp_.begin(), stamp_.end(), 0);
            now_ = 1;
        }
    }
    uint32_t dist(uint32_t v) const { return stamp_[v] == now_ ? dist_[v] : kInf; }
    uint32_t parent(uint32_t v) const { return stamp_[v] == now_ ? parent_[v] : kNone; }
    void set(uint32_t v, uint32_t d, uint32_t p) {
        stamp_[v] = now_;
        dist_[v] = d;
        parent_[v] = p;
    }

private:
    std::vector<uint32_t> dist_, parent_, stamp_;
    uint32_t now_ = 0;
};

// Plain Dijkstra. parent() is the edge id used to reach a node (forward
// search) or the reverse-edge slot (backward search).
class Dijkstra {
public:
    explicit Dijkstra(const Graph& g) : g_(g), lab_(g.n) {}

    // stops when t is settled; t = kNone runs to exhaustion
    uint32_t run(uint32_t s, uint32_t t, QueryStats* st = nullptr, bool backward = false) {
        lab_.reset();
        heap_.clear();
        lab_.set(s, 0, kNone);
        heap_.push(0, s);
        while (!heap_.empty()) {
            auto [d, u] = heap_.pop();
            if (d != lab_.dist(u)) continue;
            if (st) st->settled++;
            if (u == t) return d;
            if (!backward) {
                for (uint32_t e = g_.first_out[u]; e < g_.first_out[u + 1]; ++e)
                    relax(g_.head[e], d + g_.weight[e], e);
            } else {
                for (uint32_t k = g_.rfirst_in[u]; k < g_.rfirst_in[u + 1]; ++k)
                    relax(g_.rtail[k], d + g_.weight[g_.redge[k]], g_.redge[k]);
            }
        }
        return t == kNone ? 0 : kInf;
    }
    uint32_t dist(uint32_t v) const { return lab_.dist(v); }
    uint32_t parent_edge(uint32_t v) const { return lab_.parent(v); }

    // edge ids from s to t after a forward run
    std::vector<uint32_t> edge_path(uint32_t t) const {
        std::vector<uint32_t> p;
        for (uint32_t v = t, e; (e = lab_.parent(v)) != kNone;) {
            p.push_back(e);
            v = tail_of(e);
        }
        std::reverse(p.begin(), p.end());
        return p;
    }

private:
    uint32_t tail_of(uint32_t e) const {
        return uint32_t(std::upper_bound(g_.first_out.begin(), g_.first_out.end(), e) -
                        g_.first_out.begin()) - 1;
    }
    void relax(uint32_t v, uint32_t d, uint32_t e) {
        if (d < lab_.dist(v)) {
            lab_.set(v, d, e);
            heap_.push(d, v);
        }
    }
    const Graph& g_;
    Labels lab_;
    LazyHeap heap_;
};

// Bidirectional Dijkstra: grow balls from both ends, stop when the two
// smallest keys add up to at least the best path seen so far.
class BidirDijkstra {
public:
    explicit BidirDijkstra(const Graph& g) : g_(g), f_(g.n), b_(g.n) {}

    uint32_t run(uint32_t s, uint32_t t, QueryStats* st = nullptr) {
        f_.reset(); b_.reset();
        hf_.clear(); hb_.clear();
        f_.set(s, 0, kNone); hf_.push(0, s);
        b_.set(t, 0, kNone); hb_.push(0, t);
        uint32_t best = s == t ? 0 : kInf;
        while (!hf_.empty() || !hb_.empty()) {
            uint64_t kf = hf_.empty() ? kInf : hf_.top_key();
            uint64_t kb = hb_.empty() ? kInf : hb_.top_key();
            if (kf + kb >= best) break;
            bool fwd = kf <= kb;
            auto& h = fwd ? hf_ : hb_;
            auto& me = fwd ? f_ : b_;
            auto& other = fwd ? b_ : f_;
            auto [d, u] = h.pop();
            if (d != me.dist(u)) continue;
            if (st) st->settled++;
            uint32_t lo = fwd ? g_.first_out[u] : g_.rfirst_in[u];
            uint32_t hi = fwd ? g_.first_out[u + 1] : g_.rfirst_in[u + 1];
            for (uint32_t k = lo; k < hi; ++k) {
                uint32_t v = fwd ? g_.head[k] : g_.rtail[k];
                uint32_t nd = d + g_.weight[fwd ? k : g_.redge[k]];
                if (nd < me.dist(v)) {
                    me.set(v, nd, k);
                    h.push(nd, v);
                }
                uint32_t od = other.dist(v);
                if (od != kInf && uint64_t(nd) + od < best) best = nd + od;
            }
        }
        return best;
    }

private:
    const Graph& g_;
    Labels f_, b_;
    LazyHeap hf_, hb_;
};

// A* with a potential pi(v) that never overestimates the remaining
// distance and is consistent: pi(u) <= w(u, v) + pi(v) for every edge.
// Then reduced costs are nonnegative and the search is exact.
template <class Potential>
class AStar {
public:
    AStar(const Graph& g, Potential& pi) : g_(g), pi_(pi), lab_(g.n) {}

    uint32_t run(uint32_t s, uint32_t t, QueryStats* st = nullptr) {
        lab_.reset();
        heap_.clear();
        pi_.set_target(t);
        lab_.set(s, 0, kNone);
        heap_.push(pi_(s), s);
        while (!heap_.empty()) {
            auto [k, u] = heap_.pop();
            uint32_t d = lab_.dist(u);
            if (k != d + pi_(u)) continue;
            if (st) st->settled++;
            if (u == t) return d;
            for (uint32_t e = g_.first_out[u]; e < g_.first_out[u + 1]; ++e) {
                uint32_t v = g_.head[e], nd = d + g_.weight[e];
                if (nd < lab_.dist(v)) {
                    lab_.set(v, nd, e);
                    heap_.push(nd + pi_(v), v);
                }
            }
        }
        return kInf;
    }

private:
    const Graph& g_;
    Potential& pi_;
    Labels lab_;
    LazyHeap heap_;
};

// Straight-line potential: distance to the target at the top speed of the
// whole graph. Edge weights are rounded up and edge lengths follow the road,
// so floor(haversine * 10 / vmax) is consistent.
class GeoPotential {
public:
    explicit GeoPotential(const Graph& g) : g_(g), cache_(g.n, 0), stamp_(g.n, 0) {
        inv_ = 10.0 / g.max_speed_mps;
    }
    void set_target(uint32_t t) {
        t_ = t;
        tlat_ = g_.lat[t] * 1e-7;
        tlon_ = g_.lon[t] * 1e-7;
        if (++now_ == 0) { std::fill(stamp_.begin(), stamp_.end(), 0); now_ = 1; }
    }
    uint32_t operator()(uint32_t v) {
        if (stamp_[v] == now_) return cache_[v];
        double d = geo::haversine(g_.lat[v] * 1e-7, g_.lon[v] * 1e-7, tlat_, tlon_);
        uint32_t p = uint32_t(std::floor(d * inv_ * (1 - 1e-12)));
        stamp_[v] = now_;
        cache_[v] = p;
        return p;
    }

private:
    const Graph& g_;
    std::vector<uint32_t> cache_, stamp_;
    uint32_t now_ = 0, t_ = 0;
    double tlat_ = 0, tlon_ = 0, inv_ = 0;
};

// ALT: A* with landmarks and the triangle inequality. For a landmark L,
// d(v, t) >= d(L, t) - d(L, v) and d(v, t) >= d(v, L) - d(t, L).
class Landmarks {
public:
    Landmarks(const Graph& g, int k, uint32_t seed = 1) : g_(g), k_(k) {
        Dijkstra dj(g);
        std::mt19937 rng(seed);
        uint32_t start = std::uniform_int_distribution<uint32_t>(0, g.n - 1)(rng);
        // farthest selection: each new landmark is the node farthest from
        // the ones chosen so far
        std::vector<uint32_t> mind(g.n, kInf);
        dj.run(start, kNone);
        uint32_t next = farthest(dj, mind, true);
        for (int i = 0; i < k; ++i) {
            ids_.push_back(next);
            from_.emplace_back(g.n);
            to_.emplace_back(g.n);
            dj.run(next, kNone, nullptr, false);
            for (uint32_t v = 0; v < g.n; ++v) from_.back()[v] = dj.dist(v);
            next = farthest(dj, mind, false);
            dj.run(ids_.back(), kNone, nullptr, true);
            for (uint32_t v = 0; v < g.n; ++v) to_.back()[v] = dj.dist(v);
        }
        cache_.assign(g.n, 0);
        stamp_.assign(g.n, 0);
    }
    const std::vector<uint32_t>& ids() const { return ids_; }

    void set_target(uint32_t t) {
        t_ = t;
        if (++now_ == 0) { std::fill(stamp_.begin(), stamp_.end(), 0); now_ = 1; }
    }
    uint32_t operator()(uint32_t v) {
        if (stamp_[v] == now_) return cache_[v];
        int64_t best = 0;
        for (int i = 0; i < k_; ++i) {
            best = std::max(best, int64_t(from_[i][t_]) - from_[i][v]);
            best = std::max(best, int64_t(to_[i][v]) - to_[i][t_]);
        }
        stamp_[v] = now_;
        cache_[v] = uint32_t(best);
        return uint32_t(best);
    }

private:
    uint32_t farthest(const Dijkstra& dj, std::vector<uint32_t>& mind, bool first) {
        uint32_t arg = 0, val = 0;
        for (uint32_t v = 0; v < g_.n; ++v) {
            uint32_t d = dj.dist(v);
            mind[v] = first ? d : std::min(mind[v], d);
            if (mind[v] != kInf && mind[v] > val) { val = mind[v]; arg = v; }
        }
        return arg;
    }
    const Graph& g_;
    int k_;
    std::vector<uint32_t> ids_;
    std::vector<std::vector<uint32_t>> from_, to_;
    std::vector<uint32_t> cache_, stamp_;
    uint32_t t_ = 0, now_ = 0;
};
