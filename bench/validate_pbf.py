"""Check the hand-written PBF reader against pyosmium (libosmium).

Both sides count nodes, ways, highway ways and their node references, and
sum ids and coordinates modulo 2^64. Every line must match.

usage: python bench/validate_pbf.py data/telangana.osm.pbf ./router
"""

import subprocess
import sys

import osmium

M = 1 << 64


class Count(osmium.SimpleHandler):
    def __init__(self):
        super().__init__()
        self.v = dict(nodes=0, ways=0, highway_ways=0, highway_refs=0,
                      node_id_sum=0, coord_sum=0, highway_ref_sum=0)

    def node(self, n):
        v = self.v
        v["nodes"] += 1
        v["node_id_sum"] = (v["node_id_sum"] + n.id) % M
        lat, lon = n.location.y, n.location.x   # fixed point, 1e-7 degrees
        v["coord_sum"] = (v["coord_sum"] + 3 * lat + lon) % M

    def way(self, w):
        v = self.v
        v["ways"] += 1
        if "highway" in w.tags:
            v["highway_ways"] += 1
            v["highway_refs"] += len(w.nodes)
            v["highway_ref_sum"] = (v["highway_ref_sum"] + sum(r.ref for r in w.nodes)) % M


def main(pbf, router):
    h = Count()
    h.apply_file(pbf, locations=False)
    out = subprocess.run([router, "pbfstats", pbf], capture_output=True, text=True, check=True).stdout
    mine = {k: int(v) for k, v in (line.split() for line in out.strip().splitlines())}
    ok = True
    for k, ref in h.v.items():
        same = mine[k] == ref
        ok &= same
        print(f"{k:16s} pyosmium {ref:>22d}  router {mine[k]:>22d}  {'ok' if same else 'MISMATCH'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
