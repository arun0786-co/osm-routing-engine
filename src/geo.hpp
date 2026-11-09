#pragma once

#include <cmath>
#include <cstdint>

namespace geo {

constexpr double kEarthRadius = 6371008.8;  // mean radius, metres

inline double to_deg(int32_t fixed) { return fixed * 1e-7; }

// Great-circle distance in metres between two points given in degrees.
inline double haversine(double lat1, double lon1, double lat2, double lon2) {
    constexpr double r = M_PI / 180.0;
    double dlat = (lat2 - lat1) * r, dlon = (lon2 - lon1) * r;
    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1 * r) * std::cos(lat2 * r) * std::sin(dlon / 2) * std::sin(dlon / 2);
    return 2 * kEarthRadius * std::asin(std::min(1.0, std::sqrt(a)));
}

// Interleave the bits of two 32-bit values (Morton / Z-order code), used
// to renumber nodes so that nearby nodes sit near each other in memory.
inline uint64_t spread(uint32_t x) {
    uint64_t v = x;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFull;
    v = (v | (v << 8)) & 0x00FF00FF00FF00FFull;
    v = (v | (v << 4)) & 0x0F0F0F0F0F0F0F0Full;
    v = (v | (v << 2)) & 0x3333333333333333ull;
    v = (v | (v << 1)) & 0x5555555555555555ull;
    return v;
}

inline uint64_t morton(int32_t lat, int32_t lon) {
    uint32_t a = uint32_t(int64_t(lat) + 900000000);
    uint32_t b = uint32_t(int64_t(lon) + 1800000000);
    return spread(a) | (spread(b) << 1);
}

}  // namespace geo
