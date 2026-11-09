"""Plot query time against Dijkstra rank.

usage: python bench/plot_rank.py bench/results/telangana_rank.csv docs/img/rank.png
"""

import csv
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = {"dijkstra": "#7f7f7f", "bidirectional": "#b07aa1", "astar": "#e15759",
          "alt": "#f28e2b", "ch": "#1f4e79"}
LABELS = {"dijkstra": "Dijkstra", "bidirectional": "bidirectional Dijkstra",
          "astar": "A* (straight line)", "alt": "ALT (16 landmarks)",
          "ch": "contraction hierarchy"}


def main(src, out):
    data = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(src)):
        data[r["algorithm"]][int(r["rank"])].append(float(r["us"]))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for algo in ["dijkstra", "bidirectional", "astar", "alt", "ch"]:
        ranks = sorted(data[algo])
        x = np.log2(ranks)
        q = np.array([np.percentile(data[algo][k], [25, 50, 75]) for k in ranks])
        ax.plot(x, q[:, 1], "-o", ms=3, color=COLORS[algo], label=LABELS[algo])
        ax.fill_between(x, q[:, 0], q[:, 2], color=COLORS[algo], alpha=0.15, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("Dijkstra rank of the target (log2)")
    ax.set_ylabel("query time (microseconds)")
    ax.set_title("Telangana road network, 877k junctions: median and interquartile range")
    ax.grid(alpha=0.3, which="both")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
