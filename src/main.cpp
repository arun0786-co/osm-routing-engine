// router: build a road graph from OSM data, contract it, query it, serve it.
//
//   router build    <in.osm.pbf> <out.graph>
//   router contract <graph> <out.ch>
//   router bench    <graph> <ch> [--queries N] [--landmarks K] [--csv file]
//   router rank     <graph> <ch> [--sources N] [--csv file]
//   router table    <graph> <ch> [--size N]
//   router serve    <graph> <ch> [--port P] [--web dir]
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <thread>

#include "build.hpp"
#include "ch.hpp"
#include "graph.hpp"
#include "search.hpp"

int serve(const Graph& g, const CH& ch, int port, const std::string& web);

namespace {

using Clock = std::chrono::steady_clock;
double secs(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::map<std::string, std::string> flags(int argc, char** argv, int from) {
    std::map<std::string, std::string> f;
    for (int i = from; i < argc; ++i)
        if (std::strncmp(argv[i], "--", 2) == 0 && i + 1 < argc) {
            // two statements: in "f[argv[i] + 2] = argv[++i]" C++17 evaluates
            // the right side first, so the key would come from the value
            std::string key = argv[i] + 2;
            f[key] = argv[++i];
        }
    return f;
}

template <class T>
T flag(const std::map<std::string, std::string>& f, const char* k, T def) {
    auto it = f.find(k);
    if (it == f.end()) return def;
    if constexpr (std::is_same_v<T, std::string>) return it->second;
    else return T(std::stod(it->second));
}


int cmd_build(int argc, char** argv) {
    if (argc < 4) return 2;
    BuildStats s;
    auto t = Clock::now();
    Graph g = build_graph(argv[2], int(std::thread::hardware_concurrency()), &s);
    double total = secs(t);
    g.save(argv[3]);
    std::printf("drivable ways      %zu\n", s.ways);
    std::printf("node references    %zu (%zu distinct nodes, %zu without coordinates)\n",
                s.refs, s.osm_nodes, s.missing_coords);
    std::printf("junctions          %u, directed edges before cleanup %u\n", s.junctions, s.edges_before);
    std::printf("graph              %u nodes, %u edges (largest strongly connected component)\n", g.n, g.m());
    std::printf("geometry points    %zu\n", g.geom_lat.size());
    std::printf("fastest edge       %.0f km/h\n", g.max_speed_mps * 3.6);
    std::printf("time               ways %.2fs, nodes %.2fs, graph %.2fs, total %.2fs\n",
                s.seconds_ways, s.seconds_nodes, s.seconds_graph, total);
    return 0;
}

// Counts and checksums over the whole file, to compare against an
// independent PBF reader (bench/validate_pbf.py uses pyosmium).
int cmd_pbfstats(int argc, char** argv) {
    if (argc < 3) return 2;
    osm::PbfFile f(argv[2]);
    int th = int(std::thread::hardware_concurrency());
    std::atomic<uint64_t> nodes{0}, ways{0}, hw{0}, hw_refs{0}, id_sum{0}, coord_sum{0}, ref_sum{0};
    f.for_each_node(th, [&](size_t, int64_t id, int32_t lat, int32_t lon) {
        nodes++;
        id_sum += uint64_t(id);
        coord_sum += uint64_t(int64_t(lat)) * 3 + uint64_t(int64_t(lon));
    });
    f.for_each_way(th, [&](size_t, const osm::Way& w) {
        ways++;
        for (auto& [k, v] : w.tags)
            if (k == "highway") {
                hw++;
                hw_refs += w.refs.size();
                uint64_t s = 0;
                for (int64_t r : w.refs) s += uint64_t(r);
                ref_sum += s;
                break;
            }
    });
    std::printf("nodes %llu\nways %llu\nhighway_ways %llu\nhighway_refs %llu\n"
                "node_id_sum %llu\ncoord_sum %llu\nhighway_ref_sum %llu\n",
                (unsigned long long)nodes, (unsigned long long)ways, (unsigned long long)hw,
                (unsigned long long)hw_refs, (unsigned long long)id_sum,
                (unsigned long long)coord_sum, (unsigned long long)ref_sum);
    return 0;
}

int cmd_contract(int argc, char** argv) {
    if (argc < 4) return 2;
    auto f = flags(argc, argv, 4);
    Graph g = Graph::load(argv[2]);
    CHParams p;
    p.witness_settle = flag(f, "settle", p.witness_settle);
    p.witness_settle_sim = flag(f, "settle-sim", p.witness_settle_sim);
    CHStats st;
    CH ch = CH::build(g, p, &st);
    ch.save(argv[3]);
    std::printf("nodes %u, road edges %u, shortcuts %llu (%.2f per edge)\n", g.n, g.m(),
                (unsigned long long)st.shortcuts, double(st.shortcuts) / g.m());
    std::printf("upward edges %zu, downward edges %zu, max level %u\n", ch.up.size(), ch.dn.size(),
                st.max_level);
    std::printf("contraction time %.1fs\n", st.seconds);
    return 0;
}

struct Algo {
    std::string name;
    std::function<uint32_t(uint32_t, uint32_t, QueryStats*)> run;
};

int cmd_bench(int argc, char** argv) {
    if (argc < 4) return 2;
    auto f = flags(argc, argv, 4);
    int queries = flag(f, "queries", 1000);
    int k = flag(f, "landmarks", 16);
    std::string csv = flag<std::string>(f, "csv", "");
    Graph g = Graph::load(argv[2]);
    CH ch = CH::load(argv[3]);
    std::printf("graph: %u nodes, %u edges; CH: %zu up + %zu down edges\n", g.n, g.m(), ch.up.size(),
                ch.dn.size());

    auto t = Clock::now();
    Landmarks lm(g, k);
    std::printf("ALT preprocessing: %d landmarks in %.2fs\n", k, secs(t));

    Dijkstra dj(g);
    BidirDijkstra bd(g);
    GeoPotential gp(g);
    AStar<GeoPotential> astar(g, gp);
    AStar<Landmarks> alt(g, lm);
    CHQuery chq(ch, g);
    std::vector<Algo> algos = {
        {"dijkstra", [&](uint32_t s, uint32_t t, QueryStats* q) { return dj.run(s, t, q); }},
        {"bidirectional", [&](uint32_t s, uint32_t t, QueryStats* q) { return bd.run(s, t, q); }},
        {"astar", [&](uint32_t s, uint32_t t, QueryStats* q) { return astar.run(s, t, q); }},
        {"alt", [&](uint32_t s, uint32_t t, QueryStats* q) { return alt.run(s, t, q); }},
        {"ch", [&](uint32_t s, uint32_t t, QueryStats* q) { return chq.run(s, t, q); }},
    };

    std::mt19937 rng(42);
    std::uniform_int_distribution<uint32_t> pick(0, g.n - 1);
    std::vector<std::pair<uint32_t, uint32_t>> pairs(queries);
    for (auto& p : pairs) p = {pick(rng), pick(rng)};

    std::vector<uint32_t> ref(queries);
    std::vector<std::vector<double>> times(algos.size());
    std::vector<uint64_t> settled(algos.size(), 0);
    std::vector<int> wrong(algos.size(), 0);
    for (size_t a = 0; a < algos.size(); ++a) {
        for (int q = 0; q < queries; ++q) {
            QueryStats st;
            auto t0 = Clock::now();
            uint32_t d = algos[a].run(pairs[q].first, pairs[q].second, &st);
            times[a].push_back(secs(t0) * 1e6);
            settled[a] += st.settled;
            if (a == 0) ref[q] = d;
            else wrong[a] += d != ref[q];
        }
    }
    // CH paths unpack into real road paths of the same length
    int bad_paths = 0;
    for (int q = 0; q < std::min(queries, 1000); ++q) {
        chq.run(pairs[q].first, pairs[q].second);
        uint64_t sum = 0;
        for (uint32_t e : chq.edge_path()) sum += g.weight[e];
        bad_paths += sum != ref[q];
    }

    std::printf("\n%d random queries\n", queries);
    std::printf("%-14s %12s %12s %12s %10s %8s\n", "algorithm", "mean us", "median us", "settled", "speedup", "wrong");
    double base = 0;
    std::ofstream out;
    if (!csv.empty()) {
        out.open(csv);
        out << "algorithm,mean_us,median_us,p99_us,mean_settled,wrong\n";
    }
    for (size_t a = 0; a < algos.size(); ++a) {
        auto v = times[a];
        std::sort(v.begin(), v.end());
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        if (a == 0) base = mean;
        std::printf("%-14s %12.1f %12.1f %12.0f %9.0fx %8d\n", algos[a].name.c_str(), mean,
                    v[v.size() / 2], double(settled[a]) / queries, base / mean, wrong[a]);
        if (out) out << algos[a].name << "," << mean << "," << v[v.size() / 2] << ","
                     << v[size_t(v.size() * 0.99)] << "," << double(settled[a]) / queries << ","
                     << wrong[a] << "\n";
    }
    std::printf("CH paths checked: %d, wrong length or broken: %d\n", std::min(queries, 1000), bad_paths);
    return 0;
}

// Dijkstra rank experiment (Sanders and Schultes): from a random source,
// the target of rank 2^r is the 2^r-th node Dijkstra settles. Plotting
// query time against rank separates local from long-distance queries.
int cmd_rank(int argc, char** argv) {
    if (argc < 4) return 2;
    auto f = flags(argc, argv, 4);
    int sources = flag(f, "sources", 200);
    std::string csv = flag<std::string>(f, "csv", "rank.csv");
    Graph g = Graph::load(argv[2]);
    CH ch = CH::load(argv[3]);
    Landmarks lm(g, 16);
    Dijkstra dj(g);
    BidirDijkstra bd(g);
    GeoPotential gp(g);
    AStar<GeoPotential> astar(g, gp);
    AStar<Landmarks> alt(g, lm);
    CHQuery chq(ch, g);
    std::mt19937 rng(7);
    std::uniform_int_distribution<uint32_t> pick(0, g.n - 1);
    std::ofstream out(csv);
    out << "source,rank,algorithm,us,settled\n";

    // settle order from a source, by running Dijkstra and sorting labels
    std::vector<uint32_t> order(g.n);
    for (int i = 0; i < sources; ++i) {
        uint32_t s = pick(rng);
        dj.run(s, kNone);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
            return std::make_pair(dj.dist(a), a) < std::make_pair(dj.dist(b), b);
        });
        for (uint32_t r = 1; (1u << r) < g.n; ++r) {
            uint32_t t = order[1u << r];
            auto time = [&](const char* name, auto&& fn) {
                QueryStats st;
                auto t0 = Clock::now();
                fn(&st);
                out << s << "," << (1u << r) << "," << name << "," << secs(t0) * 1e6 << ","
                    << st.settled << "\n";
            };
            time("dijkstra", [&](QueryStats* q) { dj.run(s, t, q); });
            time("bidirectional", [&](QueryStats* q) { bd.run(s, t, q); });
            time("astar", [&](QueryStats* q) { astar.run(s, t, q); });
            time("alt", [&](QueryStats* q) { alt.run(s, t, q); });
            time("ch", [&](QueryStats* q) { chq.run(s, t, q); });
        }
    }
    std::printf("wrote %s\n", csv.c_str());
    return 0;
}

int cmd_table(int argc, char** argv) {
    if (argc < 4) return 2;
    auto f = flags(argc, argv, 4);
    int size = flag(f, "size", 1000);
    Graph g = Graph::load(argv[2]);
    CH ch = CH::load(argv[3]);
    CHQuery chq(ch, g);
    std::mt19937 rng(11);
    std::uniform_int_distribution<uint32_t> pick(0, g.n - 1);
    std::vector<uint32_t> pts(size);
    for (auto& p : pts) p = pick(rng);
    auto t = Clock::now();
    auto D = chq.table(pts, pts);
    double tt = secs(t);
    // spot check against one-to-one Dijkstra, and time the naive way
    Dijkstra dj(g);
    int checks = 0, bad = 0;
    t = Clock::now();
    for (int i = 0; i < std::min(size, 20); ++i) {
        dj.run(pts[i], kNone);
        for (int j = 0; j < size; ++j, ++checks) bad += D[size_t(i) * size + j] != dj.dist(pts[j]);
    }
    double per_row = secs(t) / std::min(size, 20);
    t = Clock::now();
    for (int i = 0; i < std::min(size, 200); ++i)
        for (int j = 0; j < 5; ++j) chq.run(pts[i], pts[(i * 7 + j) % size]);
    double per_pair = secs(t) / (std::min(size, 200) * 5);
    std::printf("%dx%d table with CH buckets: %.3fs\n", size, size, tt);
    std::printf("same with one Dijkstra per row:      %.1fs (estimated from %d rows)\n",
                per_row * size, std::min(size, 20));
    std::printf("same with one CH query per pair:     %.1fs (estimated)\n", per_pair * size * size);
    std::printf("checked %d entries against Dijkstra, %d wrong\n", checks, bad);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: router build|contract|bench|rank|table|serve ...\n");
        return 2;
    }
    std::string cmd = argv[1];
    int rc = 2;
    try {
        if (cmd == "build") rc = cmd_build(argc, argv);
        else if (cmd == "pbfstats") rc = cmd_pbfstats(argc, argv);
        else if (cmd == "contract") rc = cmd_contract(argc, argv);
        else if (cmd == "bench") rc = cmd_bench(argc, argv);
        else if (cmd == "rank") rc = cmd_rank(argc, argv);
        else if (cmd == "table") rc = cmd_table(argc, argv);
        else if (cmd == "serve" && argc >= 4) {
            auto f = flags(argc, argv, 4);
            Graph g = Graph::load(argv[2]);
            CH ch = CH::load(argv[3]);
            rc = serve(g, ch, flag(f, "port", 8080), flag<std::string>(f, "web", "web"));
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    if (rc == 2) std::fprintf(stderr, "bad arguments for %s\n", cmd.c_str());
    return rc;
}
