#include "build.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string_view>

#include "geo.hpp"

namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}

std::string_view tag(const osm::Tags& tags, std::string_view key) {
    for (auto& [k, v] : tags)
        if (k == key) return v;
    return {};
}

double default_speed(std::string_view hw) {
    struct S { std::string_view name; double kmh; };
    static const S table[] = {
        {"motorway", 90}, {"motorway_link", 45}, {"trunk", 75}, {"trunk_link", 40},
        {"primary", 60}, {"primary_link", 30}, {"secondary", 50}, {"secondary_link", 25},
        {"tertiary", 40}, {"tertiary_link", 20}, {"unclassified", 30}, {"residential", 25},
        {"living_street", 10}, {"service", 15}, {"road", 25},
    };
    for (auto& s : table)
        if (s.name == hw) return s.kmh;
    return 0;
}

// "60", "60 km/h", "40 mph"; anything else (e.g. "IN:urban") is ignored
double parse_maxspeed(std::string_view s) {
    if (s.empty() || !(s[0] >= '0' && s[0] <= '9')) return 0;
    std::string str(s);
    char* end = nullptr;
    double v = std::strtod(str.c_str(), &end);
    std::string_view rest(end);
    while (!rest.empty() && rest[0] == ' ') rest.remove_prefix(1);
    if (rest == "mph") v *= 1.609344;
    else if (!rest.empty() && rest != "km/h" && rest != "kmh") return 0;
    return (v >= 5 && v <= 130) ? v : 0;
}

bool is_one_of(std::string_view v, std::initializer_list<std::string_view> xs) {
    for (auto x : xs)
        if (v == x) return true;
    return false;
}

struct RoutableWay {
    std::vector<int64_t> refs;
    double speed_kmh;
    bool fwd, bwd;
};

}  // namespace

RoadClass classify(const osm::Tags& tags) {
    RoadClass rc;
    auto hw = tag(tags, "highway");
    double speed = default_speed(hw);
    if (speed == 0) return rc;
    if (tag(tags, "area") == "yes") return rc;

    auto mv = tag(tags, "motor_vehicle"), mc = tag(tags, "motorcar");
    auto allowed = [](std::string_view v) { return is_one_of(v, {"yes", "designated", "permissive", "destination"}); };
    auto denied = [](std::string_view v) { return is_one_of(v, {"no", "private", "agricultural", "forestry", "delivery"}); };
    if (denied(mc) || denied(mv)) return rc;
    if (denied(tag(tags, "access")) && !allowed(mc) && !allowed(mv)) return rc;

    if (hw == "service" && is_one_of(tag(tags, "service"), {"parking_aisle", "drive-through"}))
        speed = 10;
    if (double ms = parse_maxspeed(tag(tags, "maxspeed")); ms > 0) speed = ms;

    rc.routable = true;
    rc.speed_kmh = speed;
    auto ow = tag(tags, "oneway");
    auto jn = tag(tags, "junction");
    if (is_one_of(ow, {"yes", "true", "1"})) rc.backward = false;
    else if (is_one_of(ow, {"-1", "reverse"})) rc.forward = false;
    else if (is_one_of(ow, {"no", "false", "0"})) {}
    else if (hw == "motorway" || hw == "motorway_link" || is_one_of(jn, {"roundabout", "circular"}))
        rc.backward = false;
    return rc;
}

Graph build_graph(const std::string& path, int threads, BuildStats* st) {
    BuildStats local;
    BuildStats& s = st ? *st : local;
    osm::PbfFile file(path);

    // pass 1: drivable ways, kept per block so the result does not depend
    // on thread scheduling
    auto t = Clock::now();
    std::vector<std::vector<RoutableWay>> per_block(file.num_blocks());
    file.for_each_way(threads, [&](size_t bi, const osm::Way& w) {
        if (w.refs.size() < 2) return;
        RoadClass rc = classify(w.tags);
        if (!rc.routable || (!rc.forward && !rc.backward)) return;
        per_block[bi].push_back({w.refs, rc.speed_kmh, rc.forward, rc.backward});
    });
    std::vector<RoutableWay> ways;
    for (auto& b : per_block)
        for (auto& w : b) ways.push_back(std::move(w));
    per_block.clear();
    s.ways = ways.size();

    // every referenced node id, with its number of uses
    std::vector<int64_t> all;
    for (auto& w : ways) all.insert(all.end(), w.refs.begin(), w.refs.end());
    s.refs = all.size();
    std::sort(all.begin(), all.end());
    std::vector<int64_t> ids;
    std::vector<uint32_t> uses;
    for (size_t i = 0; i < all.size();) {
        size_t j = i;
        while (j < all.size() && all[j] == all[i]) ++j;
        ids.push_back(all[i]);
        uses.push_back(uint32_t(j - i));
        i = j;
    }
    all.clear();
    all.shrink_to_fit();
    auto index_of = [&](int64_t id) {
        return uint32_t(std::lower_bound(ids.begin(), ids.end(), id) - ids.begin());
    };
    // way ends are junctions too
    std::vector<uint8_t> junction(ids.size(), 0);
    for (size_t i = 0; i < ids.size(); ++i) junction[i] = uses[i] >= 2;
    for (auto& w : ways) {
        junction[index_of(w.refs.front())] = 1;
        junction[index_of(w.refs.back())] = 1;
    }
    s.seconds_ways = since(t);

    // pass 2: coordinates of the referenced nodes
    t = Clock::now();
    std::vector<int32_t> lat(ids.size(), 0), lon(ids.size(), 0);
    std::vector<uint8_t> have(ids.size(), 0);
    file.for_each_node(threads, [&](size_t, int64_t id, int32_t la, int32_t lo) {
        auto it = std::lower_bound(ids.begin(), ids.end(), id);
        if (it != ids.end() && *it == id) {
            size_t k = size_t(it - ids.begin());
            lat[k] = la;
            lon[k] = lo;
            have[k] = 1;
        }
    });
    s.osm_nodes = ids.size();
    for (auto h : have) s.missing_coords += !h;
    s.seconds_nodes = since(t);

    // junction numbering and edges between consecutive junctions
    t = Clock::now();
    std::vector<uint32_t> node_of(ids.size(), kNone);
    std::vector<int32_t> nlat, nlon;
    for (size_t i = 0; i < ids.size(); ++i)
        if (junction[i] && have[i]) {
            node_of[i] = uint32_t(nlat.size());
            nlat.push_back(lat[i]);
            nlon.push_back(lon[i]);
        }
    s.junctions = uint32_t(nlat.size());

    std::vector<EdgeIn> edges;
    for (auto& w : ways) {
        double mps = w.speed_kmh / 3.6;
        uint32_t start = kNone;
        double len = 0;
        std::vector<int32_t> glat, glon;
        size_t prev = kNone;
        for (int64_t ref : w.refs) {
            size_t k = index_of(ref);
            if (!have[k]) {  // the extract cut this way; restart at the next junction
                start = kNone;
                prev = kNone;
                continue;
            }
            if (prev != kNone && start != kNone)
                len += geo::haversine(lat[prev] * 1e-7, lon[prev] * 1e-7, lat[k] * 1e-7, lon[k] * 1e-7);
            if (node_of[k] != kNone) {
                if (start != kNone && node_of[k] != start) {
                    uint32_t wds = std::max<uint32_t>(1, uint32_t(std::ceil(len / mps * 10.0)));
                    if (w.fwd) edges.push_back({start, node_of[k], wds, float(len), glat, glon});
                    if (w.bwd) {
                        std::vector<int32_t> rl(glat.rbegin(), glat.rend()), ro(glon.rbegin(), glon.rend());
                        edges.push_back({node_of[k], start, wds, float(len), std::move(rl), std::move(ro)});
                    }
                }
                start = node_of[k];
                len = 0;
                glat.clear();
                glon.clear();
            } else if (start != kNone) {
                glat.push_back(lat[k]);
                glon.push_back(lon[k]);
            }
            prev = k;
        }
    }
    s.edges_before = uint32_t(edges.size());
    ways.clear();
    Graph g = graph_from_edges(uint32_t(nlat.size()), std::move(nlat), std::move(nlon), edges);
    s.seconds_graph = since(t);
    return g;
}
