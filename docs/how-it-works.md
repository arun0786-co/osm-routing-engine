# How it works

The pipeline has three stages. A raw OpenStreetMap extract becomes a road graph (`router build`). The graph is preprocessed into a contraction hierarchy (`router contract`). Queries then run against both (`router bench`, `router serve`). This document follows that order and ends with the bugs found on the way.

## 1. Reading OpenStreetMap data

An `.osm.pbf` file is a sequence of blobs. Each blob has three parts:
- a 4-byte big-endian length,
- a `BlobHeader` protobuf message giving the blob type and size,
- a `Blob` message whose payload is usually zlib-compressed.

The payload of an `OSMData` blob is a `PrimitiveBlock`. It holds:
- a string table, which tags refer to by index,
- groups of nodes, ways and relations,
- a coordinate granularity and offset.

There is no generated protobuf code. `ProtoReader` (`src/pbf.cpp`) decodes the wire format directly: varints, zigzag-encoded signed integers, and length-delimited fields. Most nodes come as `DenseNodes`. There, ids, latitudes and longitudes are three packed arrays, each delta-coded against the previous entry. Way node references are delta-coded the same way.

The file is memory-mapped. Blobs are independent, so they are decompressed and decoded on all cores, and each callback is told which block it came from. The builder keeps its results per block, so the output does not depend on thread scheduling. Telangana (104 MB, 13.4M nodes, 2.3M ways) is read in about 0.3 s on 12 cores.

`bench/validate_pbf.py` reads the same file with pyosmium (libosmium) and compares node, way and reference counts plus 64-bit sums of ids and coordinates. All seven numbers match.

## 2. From OSM to a road graph

**Car profile** (`classify` in `src/build.cpp`):
- A way is drivable if its `highway` tag is a road class and access is not denied for cars.
- Speed comes from the road class (motorway 90 km/h down to living street 10 km/h), or from a parseable `maxspeed` tag.
- Direction comes from `oneway` (including `-1`). Motorways, motorway links and roundabouts are one-way by default.

**Junctions.** OSM ways are polylines with many shape points. Only junctions matter for routing: nodes used by two or more ways, plus the ends of every way. The builder counts references to each node id and walks each way, emitting one edge between consecutive junctions. The shape points in between become geometry for drawing the route. For Telangana this turns 3.6M road nodes into 882k junctions.

**Weights.** Each edge's travel time is its length divided by its speed, in tenths of a second, **rounded up**. Rounding up keeps every edge at least as slow as the straight-line bound A* uses (section 4).

**Cleanup.**
- Self-loops are dropped, and of parallel edges only the fastest is kept.
- Only the largest strongly connected component is kept, so every query has an answer. Telangana keeps 876,729 of 882,421 junctions; the rest are islands such as private compounds with one-way exits.
- Nodes are renumbered along a Z-order (Morton) curve, so nodes that are close on the map are close in memory. Searches then touch fewer cache lines.

The graph is stored in compressed sparse row form: `first_out[v]..first_out[v+1]` indexes v's outgoing edges in flat `head`, `weight` and `length` arrays. A reverse copy serves backward searches.

## 3. Dijkstra and friends

Every search reuses its arrays between queries. Each node's label carries a timestamp, so "reset" just increments a counter. A query then costs time in proportion to what it touches, not to the size of the graph. The priority queue is a binary heap with lazy deletion: instead of a decrease-key, a new entry is pushed and stale entries are skipped when popped. Entries are packed into one 64-bit word (key in the high half, node in the low half), so comparisons are single integer compares.

- **Dijkstra** settles nodes in order of distance. A long query in Telangana settles about 440k nodes, half the graph, in about 24 ms.
- **Bidirectional Dijkstra** searches forward from s and backward from t on the reverse graph. It stops when the two smallest keys add up to at least the best path found so far. On a flat map, two balls of half the radius cover about half the area of one full ball. Here it settles 24% fewer nodes (334k against 441k).
- **A\*** adds a potential pi(v), a lower bound on the remaining distance, to each key. If pi is consistent (`pi(u) <= w(u,v) + pi(v)` for every edge), A* is Dijkstra on reduced costs `w(u,v) - pi(u) + pi(v) >= 0`, so it stays exact. The straight-line potential is `floor(haversine(v, t) * 10 / v_max)`. It is consistent because:
  - each edge's length follows the road, so it is at least the straight-line distance between its ends;
  - each weight is rounded up;
  - the triangle inequality holds on the sphere.

  It is weak here because v_max must be the fastest edge anywhere in the state: a 120 km/h `maxspeed` tag, while most roads are 25 to 60 km/h.
- **ALT** (A*, Landmarks, Triangle inequality; Goldberg and Harrelson 2005) precomputes exact distances from and to 16 landmarks. The landmarks are picked greedily, each as far as possible from the ones already chosen. Then `d(v,t) >= d(L,t) - d(L,v)` and `d(v,t) >= d(v,L) - d(t,L)` give a much tighter consistent potential. It settles about 9 times fewer nodes than Dijkstra.

## 4. Contraction hierarchies

The idea (Geisberger, Sanders, Schultes and Delling, 2008): order the nodes by importance, then remove them one at a time from least to most important. When node v is removed, the two-hop path u -> v -> x is shortest only through v. To preserve it, add a **shortcut** u -> x with weight `w(u,v) + w(v,x)`. Skip the shortcut if a **witness** exists: another path from u to x, avoiding v, that is no longer.

When every node has been removed, the original edges plus the shortcuts have a key property. Between any s and t there is a shortest path that first only goes **up** in importance, then only **down**. A query is therefore two Dijkstra searches that only follow upward edges: one from s on the forward graph, one from t on the backward graph. They meet at the most important node of the path. Each upward search space is tiny: about 200 settled nodes instead of 440k.

**Order** (`src/ch.cpp`, `Contractor`). The priority of a node is

    2 * (shortcuts it would add - edges it would remove) + contracted neighbours + level

The first term (edge difference) keeps the graph sparse. The second spreads contraction evenly over the map. The third keeps the hierarchy shallow. Priorities are kept in a heap with lazy updates. When a node reaches the top, its priority is recomputed; if it is no longer the smallest, it goes back in. After each contraction, the neighbours' priorities are recomputed.

**Witness search.** Before deciding the shortcuts for an incoming neighbour u, one Dijkstra runs from u in the remaining graph, skipping v. It stops at the largest `w(u,v) + w(v,x)` or after a fixed number of settled nodes: 600 when contracting, 60 when only estimating a priority. Stopping early can only miss witnesses. That adds shortcuts that were not needed, but never wrong ones, so the result is always exact and only the size of the hierarchy depends on these limits.

**Query.** Both upward searches alternate by smallest key. Each one stops when its smallest key reaches the best meeting distance found. **Stall-on-demand** prunes further: if a more important neighbour already reaches v more cheaply through a downward edge, v's label cannot lie on a shortest up-path, so its edges are not relaxed.

**Path unpacking.** Every shortcut stores the node it skipped. Shortcut u -> x via m is replaced by u -> m and m -> x. Those two edges are still in the hierarchy, stored at m, because m was contracted before both u and x. The replacement repeats recursively until only road edges remain. Their geometry gives the line drawn on the map.

**Many-to-many tables** (Knopp et al. 2007):
1. From each target, run a backward upward search. Leave `(target, distance)` in a bucket at every node it settles.
2. From each source, run a forward upward search and scan the buckets of the nodes it settles. Every shortest path meets at its top node, so the minimum over bucket entries is exact.

A 5000 x 5000 table needs 10,000 small searches instead of 25 million queries.

## 5. Bugs found on the way

1. **Flag parsing.** `f[argv[i] + 2] = argv[++i]` read the flag name from the wrong argument. Since C++17 the right-hand side of an assignment is evaluated first, so `++i` ran before the key was computed, and every flag was stored under the wrong name. Nothing failed: the default values happened to match the first runs. It showed up only when `--size 5000` printed a 1000 x 1000 table.
2. **Rounding and A\*.** Edge weights are rounded up, not to the nearest value. Otherwise an edge could be faster than the straight-line bound and A* would no longer be guaranteed exact. The tests compare A* with Dijkstra on every query.

## References

- R. Geisberger, P. Sanders, D. Schultes, D. Delling. Contraction Hierarchies: Faster and Simpler Hierarchical Routing in Road Networks. WEA 2008.
- A. V. Goldberg, C. Harrelson. Computing the Shortest Path: A* Search Meets Graph Theory. SODA 2005.
- S. Knopp, P. Sanders, D. Schultes, F. Schulz, D. Wagner. Computing Many-to-Many Shortest Paths Using Highway Hierarchies. ALENEX 2007.
- P. Sanders, D. Schultes. Highway Hierarchies Hasten Exact Shortest Path Queries. ESA 2005 (the Dijkstra-rank methodology).
- OpenStreetMap PBF format: https://wiki.openstreetmap.org/wiki/PBF_Format
