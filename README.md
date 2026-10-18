# all-paths

An exercise exploring how **paths** and **walks** behave on an urban road network. It compares how many routes exist between two nodes as the allowed route length grows. The study area is the São Paulo city portion of the Tamanduateí River basin.

- A **path** never repeats a node.
- A **walk** may revisit nodes and edges.

## What it does

`main.py` runs the full pipeline:

1. **Basin clip.** Keeps only the part of the Tamanduateí basin that lies inside São Paulo city. The city boundary is the union of the São Paulo OD zones.
2. **Study area.** Merges the OD (origin–destination survey) zones that intersect that portion. Only São Paulo zones are kept, and a zone needs at least 1% of its area inside the basin portion. The 1% rule drops slivers caused by boundary digitization mismatches.
3. **Road network.** Downloads the `drive` network from OpenStreetMap with osmnx. It keeps the largest strongly connected component as a directed graph, so one-way streets are respected and parallel edges are collapsed.
4. **Node pair.** Picks a random origin and a random destination exactly `L` edges apart (default `L = 5`).
5. **Paths.** Enumerates every simple path from origin to destination with `L` to `L + extra` edges (default `extra = 10`, so 5–15 edges).
6. **Walks.** Counts every walk in the same length range.
7. **Scatterplot.** Plots the number of paths and walks by length.
8. **Maps.** Draws one map for paths and one for walks, coloring each edge by how often the routes use it.

Lengths are measured in **number of edges**, not meters.

## Method notes

- **Walks are counted, not listed.** Walk counts grow combinatorially, so they come from sparse adjacency-matrix products:
  - `F_i = e_s·Aⁱ` counts walks from the origin with *i* edges.
  - `B_j = Aʲ·e_t` counts walks into the destination with *j* edges.
  - `F_k[t]` is the number of walks of length *k*.
  - The usage of edge `(a, b)` is `Σ_k Σ_{i<k} F_i[a]·B_{k-1-i}[b]`.

  The result is exact without enumerating anything. An edge traversed twice in one walk counts twice.
- **Paths are listed in parallel.** The search is split across the origin's successors with `joblib`. The `nx-parallel` backend is installed and enabled, but it does not implement `all_simple_paths`, so the parallelism here is explicit.

## Data

Input files live in `data-raw/`:

| File | Contents |
|---|---|
| `microbacias_tamanduatei.geojson` | Tamanduateí River micro-basins (CRS84) |
| `survey_zones.geojson` | São Paulo Metropolitan Region OD survey zones (lon/lat), with a `municipality_name` field |

## Setup

Requires Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync
```

## Usage

```bash
uv run python main.py
```

| Option | Default | Meaning |
|---|---|---|
| `--sp-length` | `5` | Shortest-path length (edges) between the two nodes |
| `--extra` | `10` | Maximum length = `sp-length + extra` |
| `--seed` | `42` | Random seed for the node pair |
| `--n-jobs` | `-1` | Parallel workers for path enumeration (`-1` = all cores) |

The first run downloads the network from the Overpass API, which takes a few minutes. The graph is then cached in `cache/`, keyed by the study-area geometry, so later runs take seconds. If Overpass answers `502 Bad Gateway`, the public server is busy: wait a bit and retry.

## Outputs

Written to `output/`:

| File | Contents |
|---|---|
| `study_area.gpkg` | Study-area polygon (EPSG:31983) |
| `counts.csv` | Number of paths and walks per length |
| `paths.csv.gz` | Every path found: length and node sequence (OSM node IDs) |
| `scatter.png` | Paths vs walks per length (log scale) |
| `map_paths.png` | Edge usage by paths (blue scale) |
| `map_walks.png` | Edge usage by walks (orange scale) |

## Tests

```bash
uv run pytest
```

The tests use a small directed grid. They check the algebraic walk counts against brute-force enumeration, and the parallel path enumeration against a direct `nx.all_simple_paths` call.
