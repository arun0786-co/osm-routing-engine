// Minimal HTTP server for the map front end. One request at a time, bound
// to 127.0.0.1. Endpoints:
//   GET /route?from=lat,lon&to=lat,lon[&algo=ch]
//   GET /compare?from=lat,lon&to=lat,lon     all algorithms on one query
//   GET /table?points=lat,lon;lat,lon;...    travel-time matrix
//   GET /, /<file>                            static files from the web dir
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "ch.hpp"
#include "geo.hpp"
#include "graph.hpp"
#include "search.hpp"

namespace {

// Uniform grid over the node coordinates for nearest-node lookups.
class NearestIndex {
public:
    explicit NearestIndex(const Graph& g, double cell_deg = 0.01) : g_(g), cell_(cell_deg) {
        for (uint32_t v = 0; v < g.n; ++v) cells_[key(cx(g.lat[v] * 1e-7), cy(g.lon[v] * 1e-7))].push_back(v);
    }
    uint32_t nearest(double lat, double lon) const {
        int x = cx(lat), y = cy(lon);
        uint32_t best = kNone;
        double bd = 1e300;
        for (int r = 0; r < 200; ++r) {
            for (int i = x - r; i <= x + r; ++i)
                for (int j = y - r; j <= y + r; ++j) {
                    if (std::max(std::abs(i - x), std::abs(j - y)) != r) continue;
                    auto it = cells_.find(key(i, j));
                    if (it == cells_.end()) continue;
                    for (uint32_t v : it->second) {
                        double d = geo::haversine(lat, lon, g_.lat[v] * 1e-7, g_.lon[v] * 1e-7);
                        if (d < bd) { bd = d; best = v; }
                    }
                }
            // everything outside ring r is at least r cells away
            if (best != kNone && bd < r * cell_ * 111000 * std::cos(lat * M_PI / 180)) break;
        }
        return best;
    }

private:
    int cx(double lat) const { return int(std::floor(lat / cell_)); }
    int cy(double lon) const { return int(std::floor(lon / cell_)); }
    static int64_t key(int x, int y) { return (int64_t(x) << 32) ^ uint32_t(y); }
    const Graph& g_;
    double cell_;
    std::unordered_map<int64_t, std::vector<uint32_t>> cells_;
};

std::string param(const std::string& query, const std::string& name) {
    size_t pos = 0;
    while (pos < query.size()) {
        size_t amp = query.find('&', pos);
        if (amp == std::string::npos) amp = query.size();
        std::string kv = query.substr(pos, amp - pos);
        size_t eq = kv.find('=');
        if (eq != std::string::npos && kv.substr(0, eq) == name) {
            std::string v = kv.substr(eq + 1), out;
            for (size_t i = 0; i < v.size(); ++i) {
                if (v[i] == '%' && i + 2 < v.size()) {
                    out += char(std::stoi(v.substr(i + 1, 2), nullptr, 16));
                    i += 2;
                } else out += v[i];
            }
            return out;
        }
        pos = amp + 1;
    }
    return "";
}

bool parse_point(const std::string& s, double& lat, double& lon) {
    return std::sscanf(s.c_str(), "%lf,%lf", &lat, &lon) == 2;
}

void reply(int fd, int code, const std::string& type, const std::string& body) {
    std::ostringstream h;
    h << "HTTP/1.1 " << code << (code == 200 ? " OK" : code == 404 ? " Not Found" : " Bad Request")
      << "\r\nContent-Type: " << type << "\r\nContent-Length: " << body.size()
      << "\r\nConnection: close\r\nCache-Control: no-store\r\n\r\n";
    std::string out = h.str() + body;
    for (size_t sent = 0; sent < out.size();) {
        ssize_t k = ::send(fd, out.data() + sent, out.size() - sent, 0);
        if (k <= 0) break;
        sent += size_t(k);
    }
}

std::string json_point(int32_t lat, int32_t lon) {
    char b[64];
    std::snprintf(b, sizeof b, "[%.6f,%.6f]", lat * 1e-7, lon * 1e-7);
    return b;
}

}  // namespace

int serve(const Graph& g, const CH& ch, int port, const std::string& web) {
    std::fprintf(stderr, "building nearest-node index and %d landmarks...\n", 16);
    NearestIndex idx(g);
    Landmarks lm(g, 16);
    Dijkstra dj(g);
    BidirDijkstra bd(g);
    GeoPotential gp(g);
    AStar<GeoPotential> astar(g, gp);
    AStar<Landmarks> alt(g, lm);
    CHQuery chq(ch, g);

    using Clock = std::chrono::steady_clock;
    auto timed = [&](const std::string& algo, uint32_t s, uint32_t t, QueryStats& st, double& us) {
        auto t0 = Clock::now();
        uint32_t d = kInf;
        if (algo == "dijkstra") d = dj.run(s, t, &st);
        else if (algo == "bidirectional") d = bd.run(s, t, &st);
        else if (algo == "astar") d = astar.run(s, t, &st);
        else if (algo == "alt") d = alt.run(s, t, &st);
        else d = chq.run(s, t, &st);
        us = std::chrono::duration<double>(Clock::now() - t0).count() * 1e6;
        return d;
    };

    int srv = ::socket(AF_INET, SOCK_STREAM, 0);
    int one = 1;
    setsockopt(srv, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(uint16_t(port));
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if (::bind(srv, reinterpret_cast<sockaddr*>(&addr), sizeof addr) < 0 || ::listen(srv, 16) < 0) {
        std::perror("bind");
        return 1;
    }
    std::fprintf(stderr, "listening on http://127.0.0.1:%d\n", port);

    while (true) {
        int fd = ::accept(srv, nullptr, nullptr);
        if (fd < 0) continue;
        char buf[16384];
        ssize_t n = ::recv(fd, buf, sizeof buf - 1, 0);
        if (n <= 0) { ::close(fd); continue; }
        buf[n] = 0;
        std::string req(buf), method, target;
        std::istringstream(req) >> method >> target;
        std::string path = target.substr(0, target.find('?'));
        std::string query = target.find('?') == std::string::npos ? "" : target.substr(target.find('?') + 1);
        try {
            if (path == "/route" || path == "/compare") {
                double la1, lo1, la2, lo2;
                if (!parse_point(param(query, "from"), la1, lo1) || !parse_point(param(query, "to"), la2, lo2)) {
                    reply(fd, 400, "application/json", "{\"error\":\"from and to must be lat,lon\"}");
                } else {
                    uint32_t s = idx.nearest(la1, lo1), t = idx.nearest(la2, lo2);
                    std::ostringstream js;
                    if (path == "/route") {
                        std::string algo = param(query, "algo");
                        if (algo.empty()) algo = "ch";
                        QueryStats st;
                        double us;
                        uint32_t d = timed(algo, s, t, st, us);
                        chq.run(s, t);  // geometry always comes from the CH path
                        double meters = 0;
                        js << "{\"algorithm\":\"" << algo << "\",\"duration_s\":" << d / 10.0
                           << ",\"query_us\":" << us << ",\"settled\":" << st.settled
                           << ",\"from\":" << json_point(g.lat[s], g.lon[s])
                           << ",\"to\":" << json_point(g.lat[t], g.lon[t]) << ",\"geometry\":[";
                        bool first = true;
                        uint32_t at = s;
                        for (uint32_t e : chq.edge_path()) {
                            meters += g.length[e];
                            auto add = [&](int32_t la, int32_t lo) {
                                js << (first ? "" : ",") << json_point(la, lo);
                                first = false;
                            };
                            add(g.lat[at], g.lon[at]);
                            for (uint32_t k = g.geom_first[e]; k < g.geom_first[e + 1]; ++k)
                                add(g.geom_lat[k], g.geom_lon[k]);
                            at = g.head[e];
                        }
                        js << (first ? "" : ",") << json_point(g.lat[at], g.lon[at]);
                        js << "],\"distance_m\":" << meters << "}";
                    } else {
                        js << "{\"results\":[";
                        const char* names[] = {"dijkstra", "bidirectional", "astar", "alt", "ch"};
                        for (int i = 0; i < 5; ++i) {
                            QueryStats st;
                            double us;
                            uint32_t d = timed(names[i], s, t, st, us);
                            js << (i ? "," : "") << "{\"algorithm\":\"" << names[i]
                               << "\",\"duration_s\":" << d / 10.0 << ",\"query_us\":" << us
                               << ",\"settled\":" << st.settled << "}";
                        }
                        js << "]}";
                    }
                    reply(fd, 200, "application/json", js.str());
                }
            } else if (path == "/table") {
                std::string pts = param(query, "points");
                std::vector<uint32_t> nodes;
                std::istringstream ps(pts);
                std::string item;
                while (std::getline(ps, item, ';')) {
                    double la, lo;
                    if (parse_point(item, la, lo)) nodes.push_back(idx.nearest(la, lo));
                }
                if (nodes.empty() || nodes.size() > 2000) {
                    reply(fd, 400, "application/json", "{\"error\":\"give 1 to 2000 points\"}");
                } else {
                    auto t0 = Clock::now();
                    auto D = chq.table(nodes, nodes);
                    double us = std::chrono::duration<double>(Clock::now() - t0).count() * 1e6;
                    std::ostringstream js;
                    js << "{\"query_us\":" << us << ",\"durations_s\":[";
                    for (size_t i = 0; i < nodes.size(); ++i) {
                        js << (i ? "," : "") << "[";
                        for (size_t j = 0; j < nodes.size(); ++j)
                            js << (j ? "," : "") << D[i * nodes.size() + j] / 10.0;
                        js << "]";
                    }
                    js << "]}";
                    reply(fd, 200, "application/json", js.str());
                }
            } else {
                std::string file = path == "/" ? "/index.html" : path;
                if (file.find("..") != std::string::npos) {
                    reply(fd, 404, "text/plain", "not found");
                } else {
                    std::ifstream in(web + file, std::ios::binary);
                    if (!in) reply(fd, 404, "text/plain", "not found");
                    else {
                        std::string body((std::istreambuf_iterator<char>(in)), {});
                        std::string type = file.size() > 5 && file.substr(file.size() - 5) == ".html"
                                               ? "text/html; charset=utf-8"
                                               : file.size() > 3 && file.substr(file.size() - 3) == ".js"
                                                     ? "text/javascript"
                                                     : "text/plain";
                        reply(fd, 200, type, body);
                    }
                }
            }
        } catch (const std::exception& e) {
            reply(fd, 400, "application/json", std::string("{\"error\":\"") + e.what() + "\"}");
        }
        ::close(fd);
    }
}
