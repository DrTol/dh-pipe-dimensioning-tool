# -*- coding: utf-8 -*-
"""
Developed on Sun Aug 31 18:53:36 2025

@author: Dr. Hakan İbrahim Tol
"""

import pandas as pd
from collections import defaultdict, deque

P = pd.read_csv("pipe_segments.csv")
req = {"id","pre_node","suc_node","length_m"}
assert req.issubset(P.columns), f"pipe_segments.csv missing {req - set(P.columns)}"

# Build adjacency and reverse maps
P = P.copy()
P["pre_node"] = P["pre_node"].astype(int)
P["suc_node"] = P["suc_node"].astype(int)
adj = defaultdict(list)      # pre_node -> list of (suc_node, edge_id)
incoming = defaultdict(list) # suc_node -> list of (pre_node, edge_id)
for r in P.itertuples(index=False):
    adj[r.pre_node].append((r.suc_node, r.id))
    incoming[r.suc_node].append((r.pre_node, r.id))

# Find single source (node in pre but not in suc)
pre = set(P["pre_node"])
suc = set(P["suc_node"])
sources = sorted(list(pre - suc))
if len(sources) != 1:
    raise SystemExit(f"[X] Need exactly ONE source; found: {sources}")
source = sources[0]

# Self-loops quick check
SL = P[P["pre_node"] == P["suc_node"]]
if len(SL):
    print("[X] Self-loops found; remove these rows:")
    print(SL)
    raise SystemExit()

# Each suc_node must have at most one predecessor (radial constraint)
dupes = P["suc_node"].value_counts()
multi_pred = dupes[dupes > 1]
if len(multi_pred):
    bad_nodes = multi_pred.index.tolist()
    print("[X] Nodes with multiple predecessors (violates radial constraint):")
    print(P[P["suc_node"].isin(bad_nodes)].sort_values("suc_node"))
    print("\nTip: keep only one incoming edge per such node (e.g., the shortest path to source).")

# DFS cycle detection with path reconstruction
visited = {}
stack = []
cycle = None
edge_by_pair = {(int(r.pre_node), int(r.suc_node)): int(r.id) for r in P.itertuples(index=False)}

def dfs(u):
    global cycle
    visited[u] = 1  # 1=visiting, 2=done
    stack.append(u)
    for v, _eid in adj[u]:
        if v not in visited:
            dfs(v)
            if cycle: return
        elif visited[v] == 1:
            # back-edge: reconstruct cycle from v -> ... -> u -> v
            i = len(stack) - 1
            while i >= 0 and stack[i] != v:
                i -= 1
            cyc_nodes = stack[i:] + [v]
            cyc_edges = []
            for a, b in zip(cyc_nodes[:-1], cyc_nodes[1:]):
                cyc_edges.append(edge_by_pair.get((a, b)))
            cycle = (cyc_nodes, cyc_edges)
            return
    stack.pop()
    visited[u] = 2

dfs(source)

if cycle:
    nodes, edges = cycle
    print("\n[X] CYCLE DETECTED")
    print("   Nodes:", " -> ".join(map(str, nodes)))
    print("   Edges:", edges)
    print("\nFix: remove ONE of the edges in the cycle (usually keep the shortest path to the source).")
else:
    print("[OK] No directed cycle found from the source.")

# Also verify all nodes are reachable and chain to source
pred_map = {}
for r in P.itertuples(index=False):
    if r.suc_node in pred_map:
        pass  # already flagged above
    pred_map[r.suc_node] = (r.pre_node, r.id, r.length_m)

def check_chain_to_source(n):
    hops = 0
    hop_cap = len(P) + len(set(P["pre_node"])) + 1
    while n != source:
        if n not in pred_map:
            raise SystemExit(f"[X] Broken path: no predecessor into node {n}.")
        n, _, _ = pred_map[n]
        hops += 1
        if hops > hop_cap:
            raise SystemExit("[X] Cycle suspected in chain-to-source check.")
    return True

for n in suc:
    check_chain_to_source(n)
