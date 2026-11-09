// Correctness tests: every algorithm must return exactly the Dijkstra
// distance, CH paths must be real paths of that length, and tables must
// match one-to-one queries.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

#include "../src/build.hpp"
#include "../src/ch.hpp"
#include "../src/geo.hpp"
#include "../src/graph.hpp"
#include "../src/pbf.hpp"
#include "../src/search.hpp"

static int failures = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

// A jittered grid of junctions about 100 m apart. Streets get random
// speeds, some are one-way and some are missing, so the graph is directed
// and irregular. Weights follow length / speed rounded up, as in the
// real builder, so the straight-line A* bound stays valid.
static Graph random_city(uint32_t side, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    uint32_t n = side * side;
    std::vector<int32_t> lat(n), lon(n);
    for (uint32_t i = 0; i < side; ++i)
        for (uint32_t j = 0; j < side; ++j) {
            lat[i * side + j] = int32_t((17.0 + i * 0.0009 + U(rng) * 0.0003) * 1e7);
            lon[i * side + j] = int32_t((78.0 + j * 0.0009 + U(rng) * 0.0003) * 1e7);
        }
    std::vector<EdgeIn> edges;
    auto road = [&](uint32_t a, uint32_t b) {
        if (U(rng) < 0.08) return;  // missing street
        double len = geo::haversine(lat[a] * 1e-7, lon[a] * 1e-7, lat[b] * 1e-7, lon[b] * 1e-7) *
                     (1.0 + 0.3 * U(rng));
        double kmh = (U(rng) < 0.1) ? 80 : 15 + 45 * U(rng);
        uint32_t w = std::max<uint32_t>(1, uint32_t(std::ceil(len / (kmh / 3.6) * 10)));
        double r = U(rng);
        if (r < 0.8 || r >= 0.9) edges.push_back({a, b, w, float(len), {}, {}});
        if (r < 0.9) edges.push_back({b, a, w, float(len), {}, {}});
    };
    for (uint32_t i = 0; i < side; ++i)
        for (uint32_t j = 0; j < side; ++j) {
            uint32_t v = i * side + j;
            if (j + 1 < side) road(v, v + 1);
            if (i + 1 < side) road(v, v + side);
            if (i + 1 < side && j + 1 < side && U(rng) < 0.1) road(v, v + side + 1);
        }
    return graph_from_edges(n, lat, lon, edges);
}

static void check_path(const Graph& g, const std::vector<uint32_t>& path, uint32_t s, uint32_t t,
                       uint32_t dist) {
    uint64_t sum = 0;
    uint32_t at = s;
    bool ok = true;
    for (uint32_t e : path) {
        uint32_t tail = uint32_t(std::upper_bound(g.first_out.begin(), g.first_out.end(), e) -
                                 g.first_out.begin()) - 1;
        if (tail != at) ok = false;
        at = g.head[e];
        sum += g.weight[e];
    }
    CHECK(ok);
    CHECK(at == t);
    CHECK(sum == dist);
}

static void test_protobuf() {
    // 300 = 0xAC 0x02; zigzag: 0 -> 0, -1 -> 1, 1 -> 2, -2 -> 3
    uint8_t buf[] = {0xAC, 0x02, 0x00, 0x01, 0x02, 0x03, 0xFE, 0xFF, 0xFF, 0xFF, 0x0F};
    osm::ProtoReader r(buf, buf + sizeof buf);
    CHECK(r.varint() == 300);
    CHECK(r.svarint() == 0);
    CHECK(r.svarint() == -1);
    CHECK(r.svarint() == 1);
    CHECK(r.svarint() == -2);
    CHECK(r.svarint() == 2147483647);
}

static void test_profile() {
    using T = osm::Tags;
    CHECK(classify(T{{"highway", "residential"}}).routable);
    CHECK(!classify(T{{"highway", "footway"}}).routable);
    CHECK(!classify(T{{"highway", "service"}, {"access", "private"}}).routable);
    CHECK(classify(T{{"highway", "service"}, {"access", "private"}, {"motor_vehicle", "yes"}}).routable);
    auto m = classify(T{{"highway", "motorway"}});
    CHECK(m.forward && !m.backward);
    auto r = classify(T{{"highway", "primary"}, {"oneway", "-1"}});
    CHECK(!r.forward && r.backward);
    auto rb = classify(T{{"highway", "tertiary"}, {"junction", "roundabout"}});
    CHECK(rb.forward && !rb.backward);
    CHECK(classify(T{{"highway", "primary"}, {"maxspeed", "40 mph"}}).speed_kmh > 64);
    CHECK(classify(T{{"highway", "primary"}, {"maxspeed", "IN:urban"}}).speed_kmh == 60);
}

static void test_algorithms(uint32_t side, uint32_t seed, int queries) {
    Graph g = random_city(side, seed);
    CHParams p;
    p.verbose = false;
    CH ch = CH::build(g, p);
    Dijkstra dj(g);
    BidirDijkstra bd(g);
    GeoPotential geo(g);
    AStar<GeoPotential> astar(g, geo);
    Landmarks lm(g, 8, seed);
    AStar<Landmarks> alt(g, lm);
    CHQuery chq(ch, g);

    std::mt19937 rng(seed + 1);
    std::uniform_int_distribution<uint32_t> pick(0, g.n - 1);
    int bad = 0;
    for (int q = 0; q < queries; ++q) {
        uint32_t s = pick(rng), t = pick(rng);
        uint32_t d = dj.run(s, t);
        uint32_t d1 = bd.run(s, t), d2 = astar.run(s, t), d3 = alt.run(s, t), d4 = chq.run(s, t);
        if (d1 != d || d2 != d || d3 != d || d4 != d) {
            if (++bad <= 5)
                std::fprintf(stderr, "s=%u t=%u dijkstra=%u bidir=%u astar=%u alt=%u ch=%u\n",
                             s, t, d, d1, d2, d3, d4);
        }
        check_path(g, chq.edge_path(), s, t, d);
        check_path(g, dj.edge_path(t), s, t, d);
    }
    CHECK(bad == 0);

    // distance table against one-to-one queries
    std::vector<uint32_t> src, dst;
    for (int i = 0; i < 25; ++i) src.push_back(pick(rng));
    for (int i = 0; i < 30; ++i) dst.push_back(pick(rng));
    auto D = chq.table(src, dst);
    int tbad = 0;
    for (size_t i = 0; i < src.size(); ++i)
        for (size_t j = 0; j < dst.size(); ++j)
            tbad += D[i * dst.size() + j] != dj.run(src[i], dst[j]);
    CHECK(tbad == 0);
    std::printf("  city %ux%u seed %u: n=%u m=%u, %d queries, %d mismatches, table errors %d\n",
                side, side, seed, g.n, g.m(), queries, bad, tbad);
}

static void test_save_load() {
    Graph g = random_city(20, 99);
    CHParams p;
    p.verbose = false;
    CH ch = CH::build(g, p);
    std::string gp = "/tmp/rt_test.graph", cp = "/tmp/rt_test.ch";
    g.save(gp);
    ch.save(cp);
    Graph g2 = Graph::load(gp);
    CH ch2 = CH::load(cp);
    CHECK(g2.n == g.n && g2.head == g.head && g2.weight == g.weight && g2.geom_lat == g.geom_lat);
    CHECK(ch2.rank == ch.rank && ch2.up_first == ch.up_first && ch2.dn.size() == ch.dn.size());
    CHQuery a(ch, g), b(ch2, g2);
    for (uint32_t s = 0; s < g.n; s += 37) CHECK(a.run(s, g.n - 1 - s % g.n) == b.run(s, g.n - 1 - s % g.n));
    std::remove(gp.c_str());
    std::remove(cp.c_str());
}

int main() {
    std::printf("protobuf\n");
    test_protobuf();
    std::printf("profile\n");
    test_profile();
    std::printf("algorithms\n");
    test_algorithms(12, 1, 2000);
    test_algorithms(40, 2, 2000);
    test_algorithms(80, 3, 1000);
    std::printf("save/load\n");
    test_save_load();
    if (failures) {
        std::printf("FAILED: %d checks\n", failures);
        return 1;
    }
    std::printf("all tests passed\n");
    return 0;
}
