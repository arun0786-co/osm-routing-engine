#include "pbf.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <atomic>
#include <stdexcept>
#include <thread>

namespace osm {

// ---------------------------------------------------------------- protobuf

uint64_t ProtoReader::varint() {
    uint64_t v = 0;
    int shift = 0;
    while (p < end) {
        uint8_t b = *p++;
        v |= uint64_t(b & 0x7f) << shift;
        if (!(b & 0x80)) return v;
        shift += 7;
        if (shift > 63) break;
    }
    throw std::runtime_error("pbf: bad varint");
}

bool ProtoReader::next(uint32_t& field, uint32_t& wire) {
    if (p >= end) return false;
    uint64_t key = varint();
    field = uint32_t(key >> 3);
    wire = uint32_t(key & 7);
    return true;
}

std::string_view ProtoReader::bytes() {
    uint64_t n = varint();
    if (n > uint64_t(end - p)) throw std::runtime_error("pbf: truncated field");
    std::string_view s(reinterpret_cast<const char*>(p), n);
    p += n;
    return s;
}

void ProtoReader::skip(uint32_t wire) {
    switch (wire) {
        case 0: varint(); break;
        case 1: p += 8; break;
        case 2: bytes(); break;
        case 5: p += 4; break;
        default: throw std::runtime_error("pbf: unsupported wire type");
    }
    if (p > end) throw std::runtime_error("pbf: truncated field");
}

// Read a repeated scalar that is normally packed but may also appear as a
// single unpacked value.
template <class F>
static void repeated(ProtoReader& r, uint32_t wire, F&& f) {
    if (wire == 2) {
        ProtoReader q(r.bytes());
        while (!q.done()) f(q);
    } else {
        f(r);
    }
}

// ---------------------------------------------------------------- file

PbfFile::PbfFile(const std::string& path) {
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) throw std::runtime_error("cannot open " + path);
    struct stat st;
    fstat(fd_, &st);
    size_ = size_t(st.st_size);
    void* m = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (m == MAP_FAILED) throw std::runtime_error("cannot mmap " + path);
    data_ = static_cast<const uint8_t*>(m);

    size_t off = 0;
    while (off + 4 <= size_) {
        uint32_t hlen = (uint32_t(data_[off]) << 24) | (uint32_t(data_[off + 1]) << 16) |
                        (uint32_t(data_[off + 2]) << 8) | uint32_t(data_[off + 3]);
        off += 4;
        if (off + hlen > size_) throw std::runtime_error("pbf: truncated header");
        ProtoReader h(data_ + off, data_ + off + hlen);
        std::string_view type;
        uint64_t datasize = 0;
        uint32_t f, w;
        while (h.next(f, w)) {
            if (f == 1 && w == 2) type = h.bytes();
            else if (f == 3 && w == 0) datasize = h.varint();
            else h.skip(w);
        }
        off += hlen;
        if (off + datasize > size_) throw std::runtime_error("pbf: truncated blob");
        if (type == "OSMData") blobs_.push_back({off, size_t(datasize)});
        else if (type != "OSMHeader") throw std::runtime_error("pbf: unknown blob type");
        off += datasize;
    }
}

PbfFile::~PbfFile() {
    if (data_) munmap(const_cast<uint8_t*>(data_), size_);
    if (fd_ >= 0) ::close(fd_);
}

void PbfFile::decode(size_t i, std::string& buf) const {
    ProtoReader r(data_ + blobs_[i].offset, data_ + blobs_[i].offset + blobs_[i].size);
    std::string_view raw, zdata;
    uint64_t raw_size = 0;
    uint32_t f, w;
    while (r.next(f, w)) {
        if (f == 1 && w == 2) raw = r.bytes();
        else if (f == 2 && w == 0) raw_size = r.varint();
        else if (f == 3 && w == 2) zdata = r.bytes();
        else if (f >= 4 && f <= 7) throw std::runtime_error("pbf: only zlib blobs are supported");
        else r.skip(w);
    }
    if (!raw.empty()) {
        buf.assign(raw);
        return;
    }
    buf.resize(raw_size);
    uLongf out = uLongf(raw_size);
    int rc = uncompress(reinterpret_cast<Bytef*>(buf.data()), &out,
                        reinterpret_cast<const Bytef*>(zdata.data()), uLong(zdata.size()));
    if (rc != Z_OK || out != raw_size) throw std::runtime_error("pbf: zlib error");
}

template <class F>
void PbfFile::parallel(int threads, F&& f) const {
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    std::exception_ptr err;
    std::atomic<bool> failed{false};
    for (int t = 0; t < std::max(1, threads); ++t) {
        pool.emplace_back([&] {
            std::string buf;
            try {
                for (size_t i; (i = next++) < blobs_.size() && !failed;) {
                    decode(i, buf);
                    f(i, std::string_view(buf));
                }
            } catch (...) {
                if (!failed.exchange(true)) err = std::current_exception();
            }
        });
    }
    for (auto& t : pool) t.join();
    if (err) std::rethrow_exception(err);
}

namespace {

struct Block {
    std::vector<std::string_view> strings;
    std::vector<std::string_view> groups;
    int64_t granularity = 100, lat_offset = 0, lon_offset = 0;

    explicit Block(std::string_view data) {
        ProtoReader r(data);
        uint32_t f, w;
        while (r.next(f, w)) {
            if (f == 1 && w == 2) {
                ProtoReader st(r.bytes());
                uint32_t sf, sw;
                while (st.next(sf, sw)) {
                    if (sf == 1 && sw == 2) strings.push_back(st.bytes());
                    else st.skip(sw);
                }
            } else if (f == 2 && w == 2) {
                groups.push_back(r.bytes());
            } else if (f == 17 && w == 0) {
                granularity = int64_t(r.varint());
            } else if (f == 19 && w == 0) {
                lat_offset = int64_t(r.varint());
            } else if (f == 20 && w == 0) {
                lon_offset = int64_t(r.varint());
            } else {
                r.skip(w);
            }
        }
    }
    // nanodegrees -> 1e-7 degrees
    int32_t lat(int64_t raw) const { return int32_t((lat_offset + granularity * raw) / 100); }
    int32_t lon(int64_t raw) const { return int32_t((lon_offset + granularity * raw) / 100); }
};

}  // namespace

void PbfFile::for_each_way(int threads,
                           const std::function<void(size_t, const Way&)>& fn) const {
    parallel(threads, [&](size_t bi, std::string_view data) {
        Block b(data);
        Way way;
        std::vector<uint32_t> keys, vals;
        for (auto g : b.groups) {
            ProtoReader gr(g);
            uint32_t f, w;
            while (gr.next(f, w)) {
                if (f != 3 || w != 2) { gr.skip(w); continue; }
                ProtoReader wr(gr.bytes());
                way.id = 0;
                way.refs.clear();
                way.tags.clear();
                keys.clear();
                vals.clear();
                uint32_t wf, ww;
                while (wr.next(wf, ww)) {
                    if (wf == 1 && ww == 0) way.id = int64_t(wr.varint());
                    else if (wf == 2) repeated(wr, ww, [&](ProtoReader& q) { keys.push_back(uint32_t(q.varint())); });
                    else if (wf == 3) repeated(wr, ww, [&](ProtoReader& q) { vals.push_back(uint32_t(q.varint())); });
                    else if (wf == 8) {
                        int64_t id = 0;
                        repeated(wr, ww, [&](ProtoReader& q) { id += q.svarint(); way.refs.push_back(id); });
                    } else wr.skip(ww);
                }
                if (keys.size() != vals.size()) throw std::runtime_error("pbf: way tag mismatch");
                for (size_t k = 0; k < keys.size(); ++k) {
                    if (keys[k] >= b.strings.size() || vals[k] >= b.strings.size())
                        throw std::runtime_error("pbf: bad string index");
                    way.tags.emplace_back(b.strings[keys[k]], b.strings[vals[k]]);
                }
                fn(bi, way);
            }
        }
    });
}

void PbfFile::for_each_node(int threads,
                            const std::function<void(size_t, int64_t, int32_t, int32_t)>& fn) const {
    parallel(threads, [&](size_t bi, std::string_view data) {
        Block b(data);
        std::vector<int64_t> ids, lats, lons;
        for (auto g : b.groups) {
            ProtoReader gr(g);
            uint32_t f, w;
            while (gr.next(f, w)) {
                if (f == 2 && w == 2) {  // DenseNodes
                    ProtoReader dr(gr.bytes());
                    ids.clear(); lats.clear(); lons.clear();
                    uint32_t df, dw;
                    while (dr.next(df, dw)) {
                        if (df == 1) {
                            int64_t v = 0;
                            repeated(dr, dw, [&](ProtoReader& q) { v += q.svarint(); ids.push_back(v); });
                        } else if (df == 8) {
                            int64_t v = 0;
                            repeated(dr, dw, [&](ProtoReader& q) { v += q.svarint(); lats.push_back(v); });
                        } else if (df == 9) {
                            int64_t v = 0;
                            repeated(dr, dw, [&](ProtoReader& q) { v += q.svarint(); lons.push_back(v); });
                        } else dr.skip(dw);
                    }
                    if (ids.size() != lats.size() || ids.size() != lons.size())
                        throw std::runtime_error("pbf: dense node arrays differ in length");
                    for (size_t k = 0; k < ids.size(); ++k) fn(bi, ids[k], b.lat(lats[k]), b.lon(lons[k]));
                } else if (f == 1 && w == 2) {  // plain Node
                    ProtoReader nr(gr.bytes());
                    int64_t id = 0, lat = 0, lon = 0;
                    uint32_t nf, nw;
                    while (nr.next(nf, nw)) {
                        if (nf == 1 && nw == 0) id = nr.svarint();
                        else if (nf == 8 && nw == 0) lat = nr.svarint();
                        else if (nf == 9 && nw == 0) lon = nr.svarint();
                        else nr.skip(nw);
                    }
                    fn(bi, id, b.lat(lat), b.lon(lon));
                } else {
                    gr.skip(w);
                }
            }
        }
    });
}

}  // namespace osm
