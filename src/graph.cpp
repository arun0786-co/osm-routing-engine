#include "graph.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <tuple>
#include <stdexcept>

#include "geo.hpp"

namespace {

// CSR adjacency from an edge list (u -> v), used for the SCC pass.
void csr(uint32_t n, const std::vector<EdgeIn>& e, bool reverse,
         std::vector<uint32_t>& first, std::vector<uint32_t>& adj) {
    first.assign(n + 1, 0);
    for (auto& x : e) first[(reverse ? x.v : x.u) + 1]++;
    for (uint32_t i = 0; i < n; ++i) first[i + 1] += first[i];
    adj.resize(e.size());
    std::vector<uint32_t> pos(first.begin(), first.end() - 1);
    for (auto& x : e) adj[pos[reverse ? x.v : x.u]++] = reverse ? x.u : x.v;
}

// Kosaraju with explicit stacks: component id of every node.
std::vector<uint32_t> scc(uint32_t n, const std::vector<EdgeIn>& e) {
    std::vector<uint32_t> f1, a1, f2, a2;
    csr(n, e, false, f1, a1);
    csr(n, e, true, f2, a2);
    std::vector<uint32_t> order;
    order.reserve(n);
    std::vector<uint8_t> seen(n, 0);
    std::vector<std::pair<uint32_t, uint32_t>> st;
    for (uint32_t s = 0; s < n; ++s) {
        if (seen[s]) continue;
        seen[s] = 1;
        st.push_back({s, f1[s]});
        while (!st.empty()) {
            auto& [v, it] = st.back();
            if (it < f1[v + 1]) {
                uint32_t w = a1[it++];
                if (!seen[w]) { seen[w] = 1; st.push_back({w, f1[w]}); }
            } else {
                order.push_back(v);
                st.pop_back();
            }
        }
    }
    std::vector<uint32_t> comp(n, kNone);
    std::vector<uint32_t> stack;
    uint32_t c = 0;
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        if (comp[*it] != kNone) continue;
        comp[*it] = c;
        stack.push_back(*it);
        while (!stack.empty()) {
            uint32_t v = stack.back();
            stack.pop_back();
            for (uint32_t k = f2[v]; k < f2[v + 1]; ++k)
                if (comp[a2[k]] == kNone) { comp[a2[k]] = c; stack.push_back(a2[k]); }
        }
        ++c;
    }
    return comp;
}

template <class T>
void write_vec(std::ofstream& f, const std::vector<T>& v) {
    uint64_t n = v.size();
    f.write(reinterpret_cast<const char*>(&n), sizeof n);
    f.write(reinterpret_cast<const char*>(v.data()), std::streamsize(n * sizeof(T)));
}

template <class T>
void read_vec(std::ifstream& f, std::vector<T>& v) {
    uint64_t n = 0;
    f.read(reinterpret_cast<char*>(&n), sizeof n);
    v.resize(n);
    f.read(reinterpret_cast<char*>(v.data()), std::streamsize(n * sizeof(T)));
    if (!f) throw std::runtime_error("graph file truncated");
}

constexpr uint64_t kMagic = 0x3148504152474452ull;  // "RDGRAPH1"

}  // namespace

Graph graph_from_edges(uint32_t n, std::vector<int32_t> lat, std::vector<int32_t> lon,
                       std::vector<EdgeIn>& edges) {
    // drop self loops, keep the fastest of parallel edges
    edges.erase(std::remove_if(edges.begin(), edges.end(),
                               [](const EdgeIn& e) { return e.u == e.v; }),
                edges.end());
    std::sort(edges.begin(), edges.end(), [](const EdgeIn& a, const EdgeIn& b) {
        return std::tie(a.u, a.v, a.w) < std::tie(b.u, b.v, b.w);
    });
    edges.erase(std::unique(edges.begin(), edges.end(),
                            [](const EdgeIn& a, const EdgeIn& b) { return a.u == b.u && a.v == b.v; }),
                edges.end());

    // keep the largest strongly connected component, so every query has an answer
    auto comp = scc(n, edges);
    std::vector<uint32_t> size;
    for (uint32_t c : comp) {
        if (c >= size.size()) size.resize(c + 1, 0);
        size[c]++;
    }
    uint32_t big = uint32_t(std::max_element(size.begin(), size.end()) - size.begin());

    // renumber the kept nodes along a Z-order curve for locality
    std::vector<uint32_t> keep;
    for (uint32_t v = 0; v < n; ++v)
        if (comp[v] == big) keep.push_back(v);
    std::sort(keep.begin(), keep.end(), [&](uint32_t a, uint32_t b) {
        return geo::morton(lat[a], lon[a]) < geo::morton(lat[b], lon[b]);
    });
    std::vector<uint32_t> id(n, kNone);
    for (uint32_t i = 0; i < keep.size(); ++i) id[keep[i]] = i;

    Graph g;
    g.n = uint32_t(keep.size());
    g.lat.resize(g.n);
    g.lon.resize(g.n);
    for (uint32_t i = 0; i < g.n; ++i) { g.lat[i] = lat[keep[i]]; g.lon[i] = lon[keep[i]]; }

    std::vector<EdgeIn*> kept;
    for (auto& e : edges)
        if (id[e.u] != kNone && id[e.v] != kNone) {
            e.u = id[e.u];
            e.v = id[e.v];
            kept.push_back(&e);
        }
    std::sort(kept.begin(), kept.end(), [](const EdgeIn* a, const EdgeIn* b) {
        return std::tie(a->u, a->v) < std::tie(b->u, b->v);
    });
    g.first_out.assign(g.n + 1, 0);
    g.geom_first.push_back(0);
    for (auto* e : kept) {
        g.first_out[e->u + 1]++;
        g.head.push_back(e->v);
        g.weight.push_back(e->w);
        g.length.push_back(e->len);
        g.geom_lat.insert(g.geom_lat.end(), e->glat.begin(), e->glat.end());
        g.geom_lon.insert(g.geom_lon.end(), e->glon.begin(), e->glon.end());
        g.geom_first.push_back(uint32_t(g.geom_lat.size()));
    }
    for (uint32_t i = 0; i < g.n; ++i) g.first_out[i + 1] += g.first_out[i];
    g.build_reverse();
    g.compute_max_speed();
    return g;
}

void Graph::build_reverse() {
    rfirst_in.assign(n + 1, 0);
    for (uint32_t e = 0; e < m(); ++e) rfirst_in[head[e] + 1]++;
    for (uint32_t i = 0; i < n; ++i) rfirst_in[i + 1] += rfirst_in[i];
    rtail.resize(m());
    redge.resize(m());
    std::vector<uint32_t> pos(rfirst_in.begin(), rfirst_in.end() - 1);
    for (uint32_t u = 0; u < n; ++u)
        for (uint32_t e = first_out[u]; e < first_out[u + 1]; ++e) {
            uint32_t k = pos[head[e]]++;
            rtail[k] = u;
            redge[k] = e;
        }
}

void Graph::compute_max_speed() {
    max_speed_mps = 0;
    for (uint32_t e = 0; e < m(); ++e)
        max_speed_mps = std::max(max_speed_mps, double(length[e]) / (weight[e] / 10.0));
}

uint32_t Graph::find_edge(uint32_t u, uint32_t v, uint32_t w) const {
    for (uint32_t e = first_out[u]; e < first_out[u + 1]; ++e)
        if (head[e] == v && weight[e] == w) return e;
    return kNone;
}

void Graph::save(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.write(reinterpret_cast<const char*>(&kMagic), sizeof kMagic);
    f.write(reinterpret_cast<const char*>(&n), sizeof n);
    write_vec(f, lat); write_vec(f, lon);
    write_vec(f, first_out); write_vec(f, head); write_vec(f, weight); write_vec(f, length);
    write_vec(f, geom_first); write_vec(f, geom_lat); write_vec(f, geom_lon);
}

Graph Graph::load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    uint64_t magic = 0;
    f.read(reinterpret_cast<char*>(&magic), sizeof magic);
    if (magic != kMagic) throw std::runtime_error(path + " is not a graph file");
    Graph g;
    f.read(reinterpret_cast<char*>(&g.n), sizeof g.n);
    read_vec(f, g.lat); read_vec(f, g.lon);
    read_vec(f, g.first_out); read_vec(f, g.head); read_vec(f, g.weight); read_vec(f, g.length);
    read_vec(f, g.geom_first); read_vec(f, g.geom_lat); read_vec(f, g.geom_lon);
    g.build_reverse();
    g.compute_max_speed();
    return g;
}
