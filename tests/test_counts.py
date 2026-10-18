from collections import Counter

import networkx as nx
import pytest

from main import count_walks, enumerate_paths


@pytest.fixture
def grid():
    G = nx.DiGraph(nx.grid_2d_graph(4, 4).to_directed())
    # Algumas vias de mão única
    G.remove_edges_from([((1, 1), (1, 2)), ((2, 2), (2, 1)), ((0, 3), (1, 3))])
    return nx.convert_node_labels_to_integers(G)


def brute_walks(G, s, t, L, Lmax):
    by_length, edge_freq = Counter(), Counter()

    def dfs(walk):
        k = len(walk) - 1
        if k >= L and walk[-1] == t:
            by_length[k] += 1
            edge_freq.update(zip(walk, walk[1:]))
        if k < Lmax:
            for n in G.successors(walk[-1]):
                dfs(walk + [n])

    dfs([s])
    return by_length, edge_freq


def test_count_walks_matches_brute_force(grid):
    s, t = 0, 15
    L = nx.shortest_path_length(grid, s, t)
    assert count_walks(grid, s, t, L, L + 4) == brute_walks(grid, s, t, L, L + 4)


def test_enumerate_paths_matches_networkx(grid):
    s, t = 0, 15
    Lmax = 9
    expected_len, expected_edges = Counter(), Counter()
    for p in nx.all_simple_paths(grid, s, t, cutoff=Lmax):
        expected_len[len(p) - 1] += 1
        expected_edges.update(zip(p, p[1:]))
    assert enumerate_paths(grid, s, t, Lmax, n_jobs=1) == (expected_len, expected_edges)


def test_adjacent_nodes(grid):
    by_length, _ = enumerate_paths(grid, 0, 1, 1, n_jobs=1)
    assert by_length == Counter({1: 1})
