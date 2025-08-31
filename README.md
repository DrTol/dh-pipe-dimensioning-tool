# District Heating Pipe Dimensioning Tool - Python
A lightweight, CSV-driven Python toolkit for hydraulic dimensioning of branched (radial) district heating (DH) networks using the Pressure Gradient / Target Pressure Loss method. The repository includes a minimal, reproducible input dataset and a QGIS project for spatial context, plus utilities for topology checks (e.g., `diagnose_cycles.py`). Results are written to a single Excel workbook with separate sheets for mains and service pipes.

> Scope: This tool targets planning- and pre-design-stage sizing for tree networks fed from a single source. It deliberately separates main vs. service connections and keeps runtime dependencies lean.

## Repository Content 
* **Core pipeline** (data preparation → network build → hydraulic sizing → export to Excel)
* **Standalone checker**: `diagnose_cycles.py` (detects graph cycles / non-radial topology)
* **QGIS project** (for map-based QA of node/edge IDs and lengths)
* **Complete example inputs** (CSV) so you can run and verify immediately

## Methodological Overview
This implementation uses a target pressure gradient (Δp/L) as the governing design variable for the supply path to the critical consumer and then checks/rightsizes all branches accordingly.

### 1. Hydraulic Headroom
The critical route length (longest hydraulic path source→leaf) is computed on the network tree. The target pressure gradient is then obtained by means of dividing the available head for distribution to the critical route length. 

### 2. Flow Aggregation on a Tree
Downstream service connections (with reference households) are summed to segment flows (space heating + DHW according to the single-row Input Data), respecting branch topology.

### 3. Diameter Selection from Catalogue
For each segment and for each service pipe, the tool iterates a friction calculation (Darcy–Weisbach with an iterative friction factor solver as implemented in the hydraulic module) and selects the smallest catalogue diameter that does not exceed the target gradient and respects velocity / regime checks.

### 4. Convergence & Checks
Iterative friction updates (configurable iteration count) stabilize the f–Re–v coupling. The algorithm reports the governing constraint (pressure gradient, velocity, or catalogue bound) for each result.

### 5. Export
Results (mains vs. services), intermediate metadata, and the original CSVs are written to one Excel workbook (separate sheets).

> **Why this method?** For branched DH in concept design, specifying a uniform permissible pressure loss per meter provides quick, consistent sizing, naturally aligned with pump head and end‐user requirements, while remaining compatible with catalogue-based diameter selection.

## How to Run

### Requirements (kept intentionally light)
Python ≥ 3.9 with (`pandas`, `numpy`, `openpyxl`)

### Run
`pipe-dimensioning-tool.py`

* `--input-dir` contains the five CSV files named exactly as above
* `--substation-type` is passed through to the hydraulic routines (for completeness in metadata)
* `--friction-iters` controls the internal friction factor update loop (stability vs. speed)

### Topology Checker
`diagnose_cycles.py`: This utility loads `pipe_segments.csv` and `nodes.csv`, constructs the directed graph, and fails fast when a cycle is present.

## Design Assumptions
* **Network type**: single-source, branched (tree) only. No loops, meshes, or multiple sources.
* **Hydraulics**: steady-state; fluid properties based on water at supplied temperatures (as implemented in the hydraulic module). Pressure drop from Darcy–Weisbach with an iterative friction factor (Colebrook-White or equivalent) and catalogue roughness.
* **Load model**: flows derived from `Space_Heating_kW` and `Domestic_Hot_Water_kW`, scaled by `ref_build` at each service point; branched aggregation toward the source.
* **Target gradient**: uniform permissible Δp/L across the network, derived from pump head and user DP via the critical route.
* **Selection logic**: choose the smallest catalogue diameter meeting Δp/L and internal velocity bounds.
* **Separation of mains vs. service pipes**: service pipes are dimensioned independently from mains with their own length and load (per service connection).

## Limitations
* **Topology**: Works only for branched networks (radial trees). No meshed grids, loops, or multiple sources/supplies.
* **Single scenario**: One global supply/return temperature pair and one global reference load row.
* **Steady-state**: No dynamic effects, storage, or transient events.
* **Single fluid**: Water properties only; glycol mixes or unusual fluids are out of scope.
* **Catalogue dependence**: Only sizes present in `pipe_catalogue.csv` can be selected.

## Citing
Tol, Hakan İbrahim. District Heating in Areas with Low Energy Houses – Detailed Analysis of District Heating Systems Based on Low Temperature Operation and Use of Renewable Energy. PhD diss. Technical University of Denmark, 2015. 
Tol, Hakan İbrahim; Svendsen, Svend. 2012. Improving the dimensioning of piping networks and network layouts in low-energy district heating systems connected to low-energy buildings: A case study in Roskilde, Denmark. Energy 38 (1): 276–290. https://doi.org/10.1016/j.energy.2011.12.002

## License
You are free to use, modify and distribute the code as long as **authorship is properly acknowledged**. Please reference this repository in derivative works.

## Acknowledgements
Above all, I give thanks to **Allah, The Creator (C.C.)**, and honor His name **Al-‘Alīm (The All-Knowing)**.

This repository is lovingly dedicated to my parents who have passed away, in remembrance of their guidance and support.

I would also like to thank **ChatGPT (by OpenAI)** for providing valuable support in updating and improving the Python implementation.
