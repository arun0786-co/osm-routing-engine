// Reader for OpenStreetMap PBF files.
//
// A PBF file is a sequence of blobs: a 4-byte big-endian length, a
// BlobHeader message, then a Blob message whose payload is usually zlib
// compressed. Each OSMData blob holds a PrimitiveBlock with a string table
// and groups of nodes, ways and relations. The protobuf wire format is
// decoded by hand here; only zlib is needed.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace osm {

using Tags = std::vector<std::pair<std::string_view, std::string_view>>;

struct Way {
    int64_t id = 0;
    std::vector<int64_t> refs;
    Tags tags;
};

// Coordinates are in units of 1e-7 degrees, the usual OSM fixed point.
struct NodeVisitor {
    virtual ~NodeVisitor() = default;
    virtual void node(int64_t id, int32_t lat, int32_t lon) = 0;
};

class PbfFile {
public:
    explicit PbfFile(const std::string& path);
    ~PbfFile();
    PbfFile(const PbfFile&) = delete;
    PbfFile& operator=(const PbfFile&) = delete;

    size_t num_blocks() const { return blobs_.size(); }
    size_t size_bytes() const { return size_; }

    // Decode all data blocks on `threads` threads. The way callback receives
    // a block index and the ways of that block; the node callback receives
    // a block index and is called once per node. Callbacks for different
    // blocks may run concurrently.
    void for_each_way(int threads,
                      const std::function<void(size_t, const Way&)>& fn) const;
    void for_each_node(int threads,
                       const std::function<void(size_t, int64_t, int32_t, int32_t)>& fn) const;

private:
    struct BlobRef { size_t offset, size; };
    void decode(size_t i, std::string& buf) const;
    template <class F> void parallel(int threads, F&& f) const;

    const uint8_t* data_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;
    std::vector<BlobRef> blobs_;
};

// Protobuf wire-format cursor, exposed for tests.
struct ProtoReader {
    const uint8_t* p;
    const uint8_t* end;

    ProtoReader(const uint8_t* b, const uint8_t* e) : p(b), end(e) {}
    explicit ProtoReader(std::string_view s)
        : p(reinterpret_cast<const uint8_t*>(s.data())), end(p + s.size()) {}

    bool done() const { return p >= end; }
    uint64_t varint();
    int64_t svarint() {
        uint64_t v = varint();
        return static_cast<int64_t>((v >> 1) ^ (~(v & 1) + 1));
    }
    // reads a field key; returns false at end of message
    bool next(uint32_t& field, uint32_t& wire);
    std::string_view bytes();
    void skip(uint32_t wire);
};

}  // namespace osm
