#include "ch.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <queue>
#include <stdexcept>

namespace {

using Adj = std::vector<std::vector<CHEdge>>;

// Dijkstra from u in the remaining graph, ignoring the node being
// contracted, cut off at a distance limit and a settle limit. Cutting it
// off early can only make it miss witnesses, which adds unneeded
// shortcuts but never wrong ones.
class Witness {
public:
    explicit Witness(uint32_t n) : lab_(n) {}
    void run(const Adj& out, uint32_t s, uint32_t skip, uint32_t limit, uint32_t max_settled) {
        lab_.reset();
        heap_.clear();
        lab_.set(s, 0, kNone);
        heap_.push(0, s);
        uint32_t settled = 0;
        while (!heap_.empty()) {
            auto [d, u] = heap_.pop();
            if (d != lab_.dist(u)) continue;
            if (d > limit || ++settled > max_settled) break;
            for (const auto& e : out[u]) {
                if (e.to == skip) continue;
                uint32_t nd = d + e.w;
                if (nd < lab_.dist(e.to)) {
                    lab_.set(e.to, nd, kNone);
                    heap_.push(nd, e.to);
                }
            }
        }
    }
    uint32_t dist(uint32_t v) const { return lab_.dist(v); }

private:
    Labels lab_;
    LazyHeap heap_;
};

class Contractor {
public:
    Contractor(const Graph& g, const CHParams& p)
        : n_(g.n), p_(p), out_(g.n), in_(g.n), done_(g.n, 0), level_(g.n, 0),
          deleted_(g.n, 0), prio_(g.n, 0), witness_(g.n) {
        for (uint32_t u = 0; u < g.n; ++u)
            for (uint32_t e = g.first_out[u]; e < g.first_out[u + 1]; ++e) {
                out_[u].push_back({g.head[e], g.weight[e], kNone});
                in_[g.head[e]].push_back({u, g.weight[e], kNone});
            }
    }

    CH run(CHStats* st) {
        auto t0 = std::chrono::steady_clock::now();
        CH ch;
        ch.n = n_;
        ch.rank.assign(n_, kNone);
        std::vector<std::vector<CHEdge>> up(n_), dn(n_);

        using Item = std::pair<int64_t, uint32_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<>> pq;
        for (uint32_t v = 0; v < n_; ++v) {
            prio_[v] = priority(v);
            pq.push({prio_[v], v});
        }
        uint64_t edges = 0;
        for (auto& o : out_) edges += o.size();
        uint32_t order = 0;
        uint32_t next_report = n_ / 20;
        while (!pq.empty()) {
            auto [p, v] = pq.top();
            pq.pop();
            if (done_[v] || p != prio_[v]) continue;
            // lazy update: the stored priority may be out of date
            int64_t np = priority(v);
            if (!pq.empty() && np > pq.top().first) {
                prio_[v] = np;
                pq.push({np, v});
                continue;
            }
            uint64_t added = contract(v, false, p_.witness_settle);
            shortcuts_ += added;
            edges += added;
            edges -= in_[v].size() + out_[v].size();
            ch.rank[v] = order++;
            done_[v] = 1;
            up[v] = std::move(out_[v]);
            dn[v] = std::move(in_[v]);
            std::vector<uint32_t> nb;
            for (auto& e : up[v]) { erase(in_[e.to], v); nb.push_back(e.to); }
            for (auto& e : dn[v]) { erase(out_[e.to], v); nb.push_back(e.to); }
            std::sort(nb.begin(), nb.end());
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            for (uint32_t x : nb) {
                deleted_[x]++;
                level_[x] = std::max(level_[x], level_[v] + 1);
                max_level_ = std::max(max_level_, level_[x]);
                prio_[x] = priority(x);
                pq.push({prio_[x], x});
            }
            out_[v].clear(); out_[v].shrink_to_fit();
            in_[v].clear(); in_[v].shrink_to_fit();
            if (p_.verbose && order >= next_report) {
                double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                std::fprintf(stderr, "  contracted %u/%u, remaining edges %llu, shortcuts %llu, %.1fs\n",
                             order, n_, (unsigned long long)edges, (unsigned long long)shortcuts_, s);
                next_report += n_ / 20;
            }
        }

        auto to_csr = [&](std::vector<std::vector<CHEdge>>& adj, std::vector<uint32_t>& first,
                          std::vector<CHEdge>& flat) {
            first.assign(n_ + 1, 0);
            for (uint32_t v = 0; v < n_; ++v) first[v + 1] = first[v] + uint32_t(adj[v].size());
            flat.reserve(first[n_]);
            for (auto& a : adj) {
                flat.insert(flat.end(), a.begin(), a.end());
                std::vector<CHEdge>().swap(a);
            }
        };
        to_csr(up, ch.up_first, ch.up);
        to_csr(dn, ch.dn_first, ch.dn);
        if (st) {
            st->shortcuts = shortcuts_;
            st->max_level = max_level_;
            st->seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        return ch;
    }

private:
    static void erase(std::vector<CHEdge>& a, uint32_t v) {
        a.erase(std::remove_if(a.begin(), a.end(), [v](const CHEdge& e) { return e.to == v; }), a.end());
    }

    // Edge difference plus terms that spread contraction evenly over the
    // graph (deleted neighbours) and keep the hierarchy shallow (level).
    int64_t priority(uint32_t v) {
        int64_t added = int64_t(contract(v, true, p_.witness_settle_sim));
        int64_t removed = int64_t(in_[v].size() + out_[v].size());
        return 2 * (added - removed) + deleted_[v] + level_[v];
    }

    uint64_t contract(uint32_t v, bool simulate, uint32_t settle) {
        uint64_t count = 0;
        for (const auto& ie : in_[v]) {
            uint32_t u = ie.to;
            uint32_t limit = 0;
            for (const auto& oe : out_[v])
                if (oe.to != u) limit = std::max(limit, ie.w + oe.w);
            if (limit == 0) continue;
            witness_.run(out_, u, v, limit, settle);
            for (const auto& oe : out_[v]) {
                uint32_t x = oe.to;
                if (x == u) continue;
                uint32_t via = ie.w + oe.w;
                if (witness_.dist(x) <= via) continue;
                ++count;
                if (!simulate) add_shortcut(u, x, via, v);
            }
        }
        return count;
    }

    void add_shortcut(uint32_t u, uint32_t x, uint32_t w, uint32_t mid) {
        for (auto& e : out_[u])
            if (e.to == x) {
                if (w < e.w) {
                    e.w = w;
                    e.mid = mid;
                    for (auto& f : in_[x])
                        if (f.to == u) { f.w = w; f.mid = mid; }
                }
                return;
            }
        out_[u].push_back({x, w, mid});
        in_[x].push_back({u, w, mid});
    }

    uint32_t n_;
    CHParams p_;
    Adj out_, in_;
    std::vector<uint8_t> done_;
    std::vector<uint32_t> level_, deleted_;
    std::vector<int64_t> prio_;
    Witness witness_;
    uint64_t shortcuts_ = 0;
    uint32_t max_level_ = 0;
};

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
    if (!f) throw std::runtime_error("ch file truncated");
}

constexpr uint64_t kMagic = 0x3148435245444f52ull;

uint32_t owner(const std::vector<uint32_t>& first, uint32_t e) {
    return uint32_t(std::upper_bound(first.begin(), first.end(), e) - first.begin()) - 1;
}

}  // namespace

CH CH::build(const Graph& g, const CHParams& p, CHStats* st) {
    Contractor c(g, p);
    return c.run(st);
}

void CH::save(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.write(reinterpret_cast<const char*>(&kMagic), sizeof kMagic);
    f.write(reinterpret_cast<const char*>(&n), sizeof n);
    write_vec(f, rank);
    write_vec(f, up_first); write_vec(f, up);
    write_vec(f, dn_first); write_vec(f, dn);
}

CH CH::load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    uint64_t magic = 0;
    f.read(reinterpret_cast<char*>(&magic), sizeof magic);
    if (magic != kMagic) throw std::runtime_error(path + " is not a CH file");
    CH ch;
    f.read(reinterpret_cast<char*>(&ch.n), sizeof ch.n);
    read_vec(f, ch.rank);
    read_vec(f, ch.up_first); read_vec(f, ch.up);
    read_vec(f, ch.dn_first); read_vec(f, ch.dn);
    return ch;
}

// ------------------------------------------------------------------ query

CHQuery::CHQuery(const CH& ch, const Graph& g) : ch_(ch), g_(g), f_(ch.n), b_(ch.n) {}

// Stall on demand: if a more important node reaches v more cheaply than
// the label being settled, v cannot be on a shortest up-path; skip it.
bool CHQuery::stalled(uint32_t v, uint32_t d, const Labels& lab, bool forward) const {
    if (forward) {
        for (uint32_t k = ch_.dn_first[v]; k < ch_.dn_first[v + 1]; ++k) {
            uint32_t du = lab.dist(ch_.dn[k].to);
            if (du != kInf && du + ch_.dn[k].w < d) return true;
        }
    } else {
        for (uint32_t k = ch_.up_first[v]; k < ch_.up_first[v + 1]; ++k) {
            uint32_t dx = lab.dist(ch_.up[k].to);
            if (dx != kInf && dx + ch_.up[k].w < d) return true;
        }
    }
    return false;
}

uint32_t CHQuery::run(uint32_t s, uint32_t t, QueryStats* st) {
    f_.reset(); b_.reset();
    hf_.clear(); hb_.clear();
    f_.set(s, 0, kNone); hf_.push(0, s);
    b_.set(t, 0, kNone); hb_.push(0, t);
    s_ = s; t_ = t;
    uint32_t best = kInf;
    meet_ = kNone;
    if (s == t) { meet_ = s; return 0; }
    while (true) {
        if (!hf_.empty() && hf_.top_key() >= best) hf_.clear();
        if (!hb_.empty() && hb_.top_key() >= best) hb_.clear();
        if (hf_.empty() && hb_.empty()) break;
        bool fwd = hb_.empty() || (!hf_.empty() && hf_.top_key() <= hb_.top_key());
        Labels& me = fwd ? f_ : b_;
        const Labels& other = fwd ? b_ : f_;
        auto [d, v] = (fwd ? hf_ : hb_).pop();
        if (d != me.dist(v)) continue;
        if (st) st->settled++;
        uint32_t od = other.dist(v);
        if (od != kInf && uint64_t(d) + od < best) { best = d + od; meet_ = v; }
        if (stalled(v, d, me, fwd)) continue;
        const auto& first = fwd ? ch_.up_first : ch_.dn_first;
        const auto& edges = fwd ? ch_.up : ch_.dn;
        for (uint32_t k = first[v]; k < first[v + 1]; ++k) {
            uint32_t x = edges[k].to, nd = d + edges[k].w;
            if (st) st->relaxed++;
            if (nd < me.dist(x)) {
                me.set(x, nd, k);
                (fwd ? hf_ : hb_).push(nd, x);
            }
        }
    }
    return best;
}

void CHQuery::unpack(uint32_t u, const CHEdge& e, std::vector<uint32_t>& out) const {
    // e is the edge u -> e.to
    if (e.mid == kNone) {
        uint32_t id = g_.find_edge(u, e.to, e.w);
        if (id == kNone) throw std::logic_error("unpack: road edge not found");
        out.push_back(id);
        return;
    }
    uint32_t m = e.mid;
    const CHEdge* a = nullptr;  // u -> m, stored at m in dn
    const CHEdge* b = nullptr;  // m -> x, stored at m in up
    for (uint32_t k = ch_.dn_first[m]; k < ch_.dn_first[m + 1]; ++k)
        if (ch_.dn[k].to == u) a = &ch_.dn[k];
    for (uint32_t k = ch_.up_first[m]; k < ch_.up_first[m + 1]; ++k)
        if (ch_.up[k].to == e.to) b = &ch_.up[k];
    if (!a || !b || a->w + b->w != e.w) throw std::logic_error("unpack: shortcut halves not found");
    unpack(u, CHEdge{m, a->w, a->mid}, out);
    unpack(m, *b, out);
}

std::vector<uint32_t> CHQuery::edge_path() const {
    std::vector<uint32_t> out;
    if (meet_ == kNone) return out;
    // source -> meeting node, along forward parents (indices into up)
    std::vector<uint32_t> fwd;
    for (uint32_t v = meet_, k; (k = f_.parent(v)) != kNone; v = owner(ch_.up_first, k)) fwd.push_back(k);
    std::reverse(fwd.begin(), fwd.end());
    for (uint32_t k : fwd) unpack(owner(ch_.up_first, k), ch_.up[k], out);
    // meeting node -> target, along backward parents (indices into dn); the
    // dn entry k stored at node y with to == v is the edge v -> y
    for (uint32_t v = meet_, k; (k = b_.parent(v)) != kNone;) {
        uint32_t y = owner(ch_.dn_first, k);
        unpack(v, CHEdge{y, ch_.dn[k].w, ch_.dn[k].mid}, out);
        v = y;
    }
    return out;
}

void CHQuery::upward(uint32_t s, Labels& lab, LazyHeap& heap, bool forward,
                     std::vector<uint32_t>* settled) {
    lab.reset();
    heap.clear();
    lab.set(s, 0, kNone);
    heap.push(0, s);
    const auto& first = forward ? ch_.up_first : ch_.dn_first;
    const auto& edges = forward ? ch_.up : ch_.dn;
    while (!heap.empty()) {
        auto [d, v] = heap.pop();
        if (d != lab.dist(v)) continue;
        if (stalled(v, d, lab, forward)) continue;
        settled->push_back(v);
        for (uint32_t k = first[v]; k < first[v + 1]; ++k) {
            uint32_t x = edges[k].to, nd = d + edges[k].w;
            if (nd < lab.dist(x)) {
                lab.set(x, nd, k);
                heap.push(nd, x);
            }
        }
    }
}

std::vector<uint32_t> CHQuery::table(const std::vector<uint32_t>& sources,
                                     const std::vector<uint32_t>& targets) {
    // backward upward search from every target, leaving (target, distance)
    // in a bucket at each node it settles; buckets are stored flat and
    // sorted by node
    struct Entry { uint32_t node, target, dist; };
    std::vector<Entry> entries;
    std::vector<uint32_t> settled;
    for (uint32_t j = 0; j < targets.size(); ++j) {
        settled.clear();
        upward(targets[j], b_, hb_, false, &settled);
        for (uint32_t v : settled) entries.push_back({v, j, b_.dist(v)});
    }
    std::sort(entries.begin(), entries.end(),
              [](const Entry& a, const Entry& b) { return a.node < b.node; });
    std::vector<uint32_t> first(ch_.n + 1, 0);
    for (auto& e : entries) first[e.node + 1]++;
    for (uint32_t v = 0; v < ch_.n; ++v) first[v + 1] += first[v];

    std::vector<uint32_t> D(sources.size() * targets.size(), kInf);
    for (uint32_t i = 0; i < sources.size(); ++i) {
        settled.clear();
        upward(sources[i], f_, hf_, true, &settled);
        uint32_t* row = &D[size_t(i) * targets.size()];
        for (uint32_t v : settled) {
            uint32_t dv = f_.dist(v);
            for (uint32_t k = first[v]; k < first[v + 1]; ++k) {
                uint32_t d = dv + entries[k].dist;
                if (d < row[entries[k].target]) row[entries[k].target] = d;
            }
        }
    }
    meet_ = kNone;
    return D;
}
