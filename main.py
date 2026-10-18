"""Paths vs walks entre dois nós da rede viária da porção paulistana da bacia do Tamanduateí."""

import argparse
import gzip
import hashlib
import tempfile
from collections import Counter
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from shapely.geometry import LineString

DATA_RAW = Path("data-raw")
CACHE = Path("cache")
OUTPUT = Path("output")
CRS_METRIC = "EPSG:31983"  # SIRGAS 2000 / UTM 23S

# Paleta (skill dataviz, modo claro)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SERIES = {"paths": "#2a78d6", "walks": "#eb6834"}
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
# Mesma luminosidade (OKLCH) da escala azul, no hue laranja da série walks
SEQ_ORANGE = ["#fbd7ca", "#f2b098", "#e68764", "#d95923", "#b44005", "#8a2e00", "#611e00"]


# 1. Porção paulistana da bacia
def sp_basin():
    basin = gpd.read_file(DATA_RAW / "microbacias_tamanduatei.geojson").to_crs(CRS_METRIC)
    zones = gpd.read_file(DATA_RAW / "survey_zones.geojson").set_crs("EPSG:4326", allow_override=True)
    zones = zones.to_crs(CRS_METRIC)
    sp_city = zones.loc[zones["municipality_name"] == "São Paulo"].union_all()
    return basin.union_all().intersection(sp_city), zones


# 2. Área de estudo = zonas OD paulistanas que intersectam a porção paulistana
# (ignora lascas de digitalização: sobreposição < min_share da área da zona)
def study_area(basin_sp, zones, min_share=0.01):
    share = zones.geometry.intersection(basin_sp).area / zones.area
    selected = zones.loc[(zones["municipality_name"] == "São Paulo") & (share >= min_share)]
    area = selected.union_all()
    OUTPUT.mkdir(exist_ok=True)
    gpd.GeoDataFrame({"name": ["study_area"]}, geometry=[area], crs=CRS_METRIC).to_file(
        OUTPUT / "study_area.gpkg"
    )
    return area, selected


# 3. Rede viária (drive)
def load_network(area):
    # Cache indexado pela geometria, para não reaproveitar a rede de outra área
    path = CACHE / f"network_{hashlib.md5(area.wkb).hexdigest()[:10]}.graphml"
    if path.exists():
        G = ox.load_graphml(path)
    else:
        area_4326 = gpd.GeoSeries([area], crs=CRS_METRIC).to_crs("EPSG:4326").iloc[0]
        G = ox.graph_from_polygon(area_4326, network_type="drive")
        CACHE.mkdir(exist_ok=True)
        ox.save_graphml(G, path)
    G = ox.truncate.largest_component(G, strongly=True)
    return nx.DiGraph(G)


# 4. Par de nós a exatamente L arestas de distância
def pick_nodes(G, L, rng):
    nodes = list(G)
    while True:
        s = nodes[rng.integers(len(nodes))]
        dist = nx.single_source_shortest_path_length(G, s, cutoff=L)
        candidates = [n for n, d in dist.items() if d == L]
        if candidates:
            return s, candidates[rng.integers(len(candidates))]


# 5. Paths simples entre s e t com comprimento <= Lmax
def _paths_from(H, s, n, t, Lmax, out_file):
    by_length, edge_freq = Counter(), Counter()
    paths = [[s, t]] if n == t else ([s, *p] for p in nx.all_simple_paths(H, n, t, cutoff=Lmax - 1))
    with gzip.open(out_file, "wt") as f:
        for p in paths:
            by_length[len(p) - 1] += 1
            edge_freq.update(zip(p, p[1:]))
            f.write(f"{len(p) - 1},{' '.join(map(str, p))}\n")
    return by_length, edge_freq


def enumerate_paths(G, s, t, Lmax, n_jobs, out_path=None):
    H = G.copy()
    H.remove_node(s)
    with tempfile.TemporaryDirectory() as tmp:
        parts = [Path(tmp) / f"part_{i}.csv.gz" for i in range(G.out_degree(s))]
        results = Parallel(n_jobs=n_jobs)(
            delayed(_paths_from)(H, s, n, t, Lmax, part) for n, part in zip(G.successors(s), parts)
        )
        if out_path is not None:
            # Membros gzip concatenados formam um gzip válido
            with open(out_path, "wb") as out:
                out.write(gzip.compress(b"length,nodes\n"))
                for part in parts:
                    out.write(part.read_bytes())
    by_length, edge_freq = Counter(), Counter()
    for bl, ef in results:
        by_length.update(bl)
        edge_freq.update(ef)
    return by_length, edge_freq


# 6. Walks entre s e t com comprimento em [L, Lmax], contados via adjacência
def count_walks(G, s, t, L, Lmax):
    nodes = list(G)
    idx = {n: i for i, n in enumerate(nodes)}
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, dtype=np.int64, format="csr")
    F = np.zeros((Lmax + 1, len(nodes)), dtype=np.int64)  # F[i, a] = walks s->a com i arestas
    B = np.zeros((Lmax + 1, len(nodes)), dtype=np.int64)  # B[j, b] = walks b->t com j arestas
    F[0, idx[s]] = 1
    B[0, idx[t]] = 1
    AT = A.T.tocsr()
    for i in range(Lmax):
        F[i + 1] = AT @ F[i]
        B[i + 1] = A @ B[i]

    by_length = Counter({k: int(F[k, idx[t]]) for k in range(L, Lmax + 1) if F[k, idx[t]]})

    rows, cols = A.nonzero()
    freq = np.zeros(len(rows), dtype=np.int64)
    for k in range(L, Lmax + 1):
        for i in range(k):
            freq += F[i, rows] * B[k - 1 - i, cols]
    edge_freq = Counter(
        {(nodes[a], nodes[b]): int(f) for a, b, f in zip(rows, cols, freq) if f}
    )
    return by_length, edge_freq


# 7. Scatterplot das duas distribuições
def _style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0)
    ax.grid(axis="y", color=GRID, linewidth=1)
    ax.set_axisbelow(True)


def plot_scatter(counts, L, out_path):
    fig, ax = plt.subplots(figsize=(8, 5), facecolor=SURFACE)
    _style_axes(ax)
    for name, color in SERIES.items():
        d = counts.loc[counts[name] > 0]
        ax.scatter(d["length"], d[name], s=64, color=color, edgecolor=SURFACE, linewidth=2,
                   zorder=3, label=name.capitalize())
        last = d.iloc[-1]
        ax.annotate(f"{name.capitalize()}: {int(last[name]):,}".replace(",", "."),
                    (last["length"], last[name]), xytext=(8, 0), textcoords="offset points",
                    va="center", color=INK_2, fontsize=9)
    ax.set_yscale("log")
    ax.set_xticks(counts["length"])
    ax.set_xlim(counts["length"].min() - 0.5, counts["length"].max() + 2.5)
    ax.set_xlabel("Comprimento (nº de arestas)", color=INK_2)
    ax.set_ylabel("Quantidade (escala log)", color=INK_2)
    ax.set_title(f"Paths vs walks por comprimento (shortest path = {L})", color=INK, loc="left")
    ax.legend(frameon=False, labelcolor=INK_2, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# 8. Mapa de uso das arestas
def _edges_gdf(G):
    rows = []
    for u, v, data in G.edges(data=True):
        geom = data.get("geometry") or LineString(
            [(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])]
        )
        rows.append({"u": u, "v": v, "geometry": geom})
    return gpd.GeoDataFrame(rows, crs="EPSG:4326").to_crs(CRS_METRIC)


def plot_map(edges, area, s, t, G, edge_freq, ramp, title, out_path):
    edges = edges.assign(freq=[edge_freq.get((u, v), 0) for u, v in zip(edges["u"], edges["v"])])
    used = edges.loc[edges["freq"] > 0].sort_values("freq")
    cmap = LinearSegmentedColormap.from_list("seq", ramp)
    norm = LogNorm(vmin=used["freq"].min(), vmax=max(used["freq"].max(), used["freq"].min() + 1))

    fig, ax = plt.subplots(figsize=(9, 9), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    gpd.GeoSeries([area], crs=CRS_METRIC).boundary.plot(ax=ax, color=MUTED, linewidth=1)
    edges.plot(ax=ax, color=GRID, linewidth=0.4)
    used.plot(ax=ax, column="freq", cmap=cmap, norm=norm, linewidth=1.6)
    st = gpd.GeoSeries(
        gpd.points_from_xy([G.nodes[n]["x"] for n in (s, t)], [G.nodes[n]["y"] for n in (s, t)]),
        crs="EPSG:4326",
    ).to_crs(CRS_METRIC)
    for pt, label in zip(st, ("Origem", "Destino")):
        ax.scatter(pt.x, pt.y, s=80, color=INK, edgecolor=SURFACE, linewidth=2, zorder=5)
        ax.annotate(label, (pt.x, pt.y), xytext=(8, 6), textcoords="offset points",
                    color=INK, fontsize=10, fontweight="bold", zorder=6)
    # Zoom na região usada, com margem
    xmin, ymin, xmax, ymax = used.total_bounds
    pad = max(xmax - xmin, ymax - ymin) * 0.15 + 300
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_axis_off()
    ax.set_title(title, color=INK, loc="left")
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Uso da aresta (escala log)", color=INK_2)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(colors=MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sp-length", type=int, default=5, help="shortest path entre os nós (arestas)")
    parser.add_argument("--extra", type=int, default=10, help="comprimento máximo = mínimo + extra")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=-1)
    args = parser.parse_args()
    nx.config.backends.parallel.n_jobs = args.n_jobs
    L, Lmax = args.sp_length, args.sp_length + args.extra

    basin_sp, zones = sp_basin()
    area, selected = study_area(basin_sp, zones)
    print(f"Porção paulistana da bacia: {basin_sp.area / 1e6:.1f} km²")
    print(f"Área de estudo: {len(selected)} zonas OD, {area.area / 1e6:.1f} km²")

    G = load_network(area)
    print(f"Rede: {G.number_of_nodes()} nós, {G.number_of_edges()} arestas")

    s, t = pick_nodes(G, L, np.random.default_rng(args.seed))
    print(f"Origem {s} -> destino {t} (shortest path = {L} arestas, máximo = {Lmax})")

    paths_len, paths_edges = enumerate_paths(G, s, t, Lmax, args.n_jobs, OUTPUT / "paths.csv.gz")
    walks_len, walks_edges = count_walks(G, s, t, L, Lmax)

    counts = pd.DataFrame({
        "length": range(L, Lmax + 1),
        "paths": [paths_len.get(k, 0) for k in range(L, Lmax + 1)],
        "walks": [walks_len.get(k, 0) for k in range(L, Lmax + 1)],
    })
    counts.to_csv(OUTPUT / "counts.csv", index=False)
    print(counts.to_string(index=False))
    print(f"Total: {counts['paths'].sum():,} paths, {counts['walks'].sum():,} walks")

    plot_scatter(counts, L, OUTPUT / "scatter.png")
    edges = _edges_gdf(G)
    plot_map(edges, area, s, t, G, paths_edges, SEQ_BLUE,
             f"Uso das arestas pelos paths ({L}–{Lmax} arestas)", OUTPUT / "map_paths.png")
    plot_map(edges, area, s, t, G, walks_edges, SEQ_ORANGE,
             f"Uso das arestas pelos walks ({L}–{Lmax} arestas)", OUTPUT / "map_walks.png")
    print(f"Saídas em {OUTPUT}/")


if __name__ == "__main__":
    main()
