#!/usr/bin/env python3
import pandas as pd
from collections import defaultdict

# Input: headerless CSV -> name, prefix, guide, suffix, target
df = pd.read_csv("guides_MPRA_dedup.csv", header=None)
df.columns = ["guide_name","prefix","guide_sequence","suffix","target"]

seqs = df["guide_sequence"].str.upper().tolist()
seq_set = set(seqs)
bases = ["A","C","G","T"]

# Build conflict graph: edge if Hamming distance 1
adj = defaultdict(set)
for s in seqs:
    for i, b in enumerate(s):
        for nb in bases:
            if nb == b:
                continue
            ns = s[:i] + nb + s[i+1:]
            if ns in seq_set:
                adj[s].add(ns)
                adj[ns].add(s)

# Greedy coloring: highest degree first
color_of = {}
for s in sorted(seqs, key=lambda x: len(adj[x]), reverse=True):
    used = {color_of[n] for n in adj[s] if n in color_of}
    c = 0
    while c in used:
        c += 1
    color_of[s] = c

n_colors = max(color_of.values()) + 1
print({
    "n_guides": len(seqs),
    "n_conflict_edges_mm1": sum(len(v) for v in adj.values()) // 2,
    "n_partitions": n_colors
})

# Write per-partition CSVs (no header, 5 columns)
for c in range(n_colors):
    keep = df[df["guide_sequence"].str.upper().map(lambda s: color_of[s]) == c]
    out = f"guides_subset_{c+1}.csv"
    keep.to_csv(out, index=False, header=False)
    print(f"wrote {out}: {len(keep)} guides")
