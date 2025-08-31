from __future__ import annotations
"""
Developed on Sun Aug 31 18:30:48 2025

@author: Dr. Hakan İbrahim Tol

"""

"""
District Heating Pipe & Service Pipe Dimensioning 
---------------------------------------------------------------------
- Input files (CSV):
    * input-data.csv (semicolon ";" separated)
      Columns:
        Space_Heating_kW;Domestic_Hot_Water_kW;Supply_Temperature_degC;Return_Temperature_sh_degC;Return_Temperature_dhw_degC;Pump_Head_Lift_bar;End-User_Differential_Pressure_bar
    * pipe_segments.csv (comma separated)
      Columns: id,pre_node,suc_node,length_m
    * nodes.csv (comma separated)
      Columns: id
    * service_pipes.csv (comma separated)
      Columns: id,node_connc,ref_build,length_m
    * Optional: pipe_catalogue.csv (comma separated)
      Columns: InnerDiameter_mm,Roughness_mm

- Output: Excel file with separate sheets (InputData, PipeList, NodeList, ServicePipes)
  where PipeList and ServicePipes include selected inner diameters.

Assumptions (adjust in code if needed):
- The network is radial (a tree). Each successor node has exactly one predecessor edge.
- Heat demands in input-data.csv are PER END-USER (reference household). "ref_build" on
  each service pipe indicates how many households connect at that service point.
- Water properties are evaluated at ~70°C (as in the legacy scripts).
- If no pipe_catalogue.csv is provided in the working directory, a default catalogue
  is used with roughness = 0.1 mm.
- Simultaneity factors reuse the legacy formulas.

Run:
  python dh_dimensioning.py  # uses files from current working directory

You can also import and call main(input_dir=Path(...)) from another script.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import math
import pandas as pd

# -----------------------------
# Hydraulic subfunctions (merged)
# -----------------------------

# Constants at ~70°C (legacy values)
G = 9.80665          # [m/s2]
RHO = 977.74         # [kg/m3]
MU = 0.0004024       # [N·s/m2]
BAR_PER_PASCAL = 1.0 / 1.0e5
PASCAL_PER_BAR = 1.0e5
BAR_PER_M_WATER = 1 / 10.1971621297792  # used in legacy PressureLoss


def reynolds(m_flow_kg_s: float, d_mm: float) -> float:
    """Reynolds number for water at ~70°C.
    m_flow_kg_s : mass flow [kg/s]
    d_mm        : inner diameter [mm]
    """
    d_m = d_mm / 1000.0
    if d_m <= 0:
        return 0.0
    # Re = 4 * m_dot / (pi * D * mu)
    return 4.0 * m_flow_kg_s / (math.pi * d_m * MU)


def f_clamond(Re: float, rel_roughness: float, iters: int = 2) -> float:
    """Efficient resolution of the Colebrook equation (Didier Clamond).

    Parameters
    ----------
    Re : float
        Reynolds number
    rel_roughness : float
        k / D (absolute roughness / diameter) [dimensionless]
    iters : int
        Number of Newton-like refinement iterations (default 2).

    Returns
    -------
    f : float
        Darcy friction factor
    """
    if Re <= 0:
        return 0.0

    # Clamond's recommended initialization 
    # X1 = K * R * log(10) / 18.574, with K = rel_roughness
    X1 = rel_roughness * Re * 0.123968186335417556
    # X2 = log( R * log(10) / 5.02 )
    X2 = math.log(Re) - 0.779397488455682028

    # Initial guess
    W = X2 - 0.2

    # Iterations
    for _ in range(max(0, iters)):
        E = (math.log(X1 + W) + W - X2) / (1.0 + X1 + W)
        W = W - (1 + X1 + W + 0.5 * E) * E * (X1 + W) / (1 + X1 + W + E * (1 + E / 3))

    # Finalized solution
    W = 1.151292546497022842 / W        # 0.5*log(10) / W
    f = W * W                           # friction factor
    return f


def pressure_loss_bar(L_m: float, d_mm: float, m_flow_kg_s: float, roughness_mm: float, iters: int = 2) -> float:
    """Pressure loss in a circular pipe segment (full flow, SI units).

    Returns pressure drop in bar.
    """
    Re = reynolds(m_flow_kg_s, d_mm)
    if Re < 2300 and Re > 0:
        f = 64.0 / Re
    elif Re <= 0:
        f = 0.0
    else:
        rel_k = (roughness_mm / max(d_mm, 1e-9))
        f = f_clamond(Re, rel_k, iters=iters)

    d_m = d_mm / 1000.0
    if d_m <= 0:
        return 0.0

    # Δp = 8 f L m_dot^2 / (pi^2 rho^2 g D^5)  (legacy formula -> returns [bar])
    dp_bar = (8 * f * L_m * m_flow_kg_s ** 2 / (math.pi ** 2 * RHO ** 2 * G * d_m ** 5)) / 10.1971621297792
    return dp_bar


# -----------------------------
# Thermo-physical helper
# -----------------------------

def h_water_kJ_per_kg(T_C: float) -> float:
    """Approx. water specific enthalpy at T (legacy linearized fit).
    Returns kJ/kg 
    """
    return 4.18683 * T_C + 0.079769  # kJ/kg


# -----------------------------
# Simultaneity factors (legacy formulas)
# -----------------------------

def sf_sh(CN: int) -> float:
    if CN <= 0:
        return 0.0
    return 0.62 + 0.38 / CN


def sf_dhw(CN: int, substation_type: int) -> float:
    if CN <= 0:
        return 0.0
    # 0: with storage; 1: without storage (legacy mapping)
    if substation_type == 0:
        A, B, C = 1.19, 1.5, 0.3
    else:
        A, B, C = 1.19, 18.0, 13.1
    peak_ref = A * 1 + B * (1 ** 0.5) + C
    return (A * CN + B * (CN ** 0.5) + C) / (peak_ref * CN)


# -----------------------------
# I/O helpers
# -----------------------------

def read_input_data(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=';')
    if df.empty:
        raise ValueError("input-data.csv is empty")
    # use first row
    return df.iloc[0]


def load_pipe_segments(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"id", "pre_node", "suc_node", "length_m"}
    if not required.issubset(df.columns):
        raise ValueError(f"pipe_segments.csv must contain columns {required}")
    return df.copy()


def load_nodes(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "id" not in df.columns:
        raise ValueError("nodes.csv must contain column 'id'")
    return df.copy()


def load_service_pipes(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"id", "node_connc", "ref_build", "length_m"}
    if not required.issubset(df.columns):
        raise ValueError(f"service_pipes.csv must contain columns {required}")
    return df.copy()


def load_pipe_catalogue(optional_path):
    if optional_path and optional_path.exists():
        try:
            cat = pd.read_csv(optional_path, sep=None, engine="python")
        except Exception:
            cat = pd.read_csv(optional_path)
        if len(cat.columns) == 1 and ";" in str(cat.columns[0]):
            cat = pd.read_csv(optional_path, sep=";")
        required = {"InnerDiameter_mm","Roughness_mm"}
        if not required.issubset(cat.columns):
            raise ValueError("pipe_catalogue.csv must contain 'InnerDiameter_mm' and 'Roughness_mm'")
        cat = cat.copy()
        cat["InnerDiameter_mm"] = pd.to_numeric(cat["InnerDiameter_mm"], errors="coerce")
        cat["Roughness_mm"]     = pd.to_numeric(cat["Roughness_mm"], errors="coerce")
        cat = cat.dropna(subset=["InnerDiameter_mm","Roughness_mm"]).sort_values("InnerDiameter_mm").reset_index(drop=True)
        return cat
    default_d = [20,25,32,40,50,65,80,100,125,150,200,250,300]
    return pd.DataFrame({"InnerDiameter_mm":default_d,"Roughness_mm":[0.1]*len(default_d)})



# -----------------------------
# Pre-processing (network adjustment)
# -----------------------------

def find_source_and_leaves(pipes: pd.DataFrame) -> Tuple[int, List[int]]:
    pre_set = set(pipes["pre_node"].astype(int))
    suc_set = set(pipes["suc_node"].astype(int))
    sources = sorted(list(pre_set - suc_set))
    if len(sources) != 1:
        raise ValueError(f"Network must have exactly one source; found {sources}.")
    source = sources[0]
    leaves = sorted([n for n in suc_set if n not in pre_set])
    if len(leaves) == 0:
        raise ValueError("No leaf nodes detected.")
    return source, leaves


def build_pred_map(pipes: pd.DataFrame) -> Dict[int, Tuple[int, int]]:
    """Map successor node -> (predecessor node, edge_index).
    Requires radial network: each suc_node appears at most once.
    """
    pred_map: Dict[int, Tuple[int, int]] = {}
    for idx, row in pipes.reset_index(drop=True).iterrows():
        suc = int(row["suc_node"])
        pre = int(row["pre_node"])
        if suc in pred_map:
            raise ValueError(f"Node {suc} has multiple predecessors; network must be radial.")
        pred_map[suc] = (pre, idx)
    return pred_map


def route_lengths_to_source(leaves, source, pred_map, lengths):
    L = []
    hop_cap = len(pred_map) + len(lengths) + 1
    for leaf in leaves:
        s = leaf
        total = 0.0
        hops = 0
        while s != source:
            if s not in pred_map:
                raise ValueError(f"Broken path: no predecessor edge into node {s}")
            pre, eidx = pred_map[s]
            total += float(lengths[eidx])
            s = pre
            hops += 1
            if hops > hop_cap:
                raise ValueError("Cycle detected while tracing path to source. Ensure the network is radial (tree).")
        L.append(total)
    return L

def cumulative_consumers_per_edge(leaves, source, pred_map, node_consumers, n_edges):
    totals = [0]*n_edges
    hop_cap = len(pred_map) + n_edges + 1
    for leaf in leaves:
        s = leaf
        cum = int(node_consumers.get(s, 0))
        hops = 0
        while s != source:
            if s not in pred_map:
                raise ValueError(f"Broken path: no predecessor edge into node {s}")
            pre, eidx = pred_map[s]
            totals[eidx] += cum
            s = pre
            cum += int(node_consumers.get(s, 0))
            hops += 1
            if hops > hop_cap:
                raise ValueError("Cycle detected while aggregating consumers. Ensure the network is radial (tree).")
    return totals

# -----------------------------
# Core dimensioning logic
# -----------------------------

def dimension_network(
    input_row: pd.Series,
    pipes: pd.DataFrame,
    nodes: pd.DataFrame,
    service: pd.DataFrame,
    catalogue: pd.DataFrame,
    substation_type: int = 1,
    friction_iters: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    """Return (pipe_results_df, service_results_df, meta).

    pipe_results_df: original pipe columns plus 'Selected_D_mm'.
    service_results_df: original service columns plus 'Selected_D_mm'.
    meta: computed values like critical_length_m and max_pressure_gradient_bar_per_m.
    """
    # Extract input
    Q_sh = float(input_row["Space_Heating_kW"])  # per end-user
    Q_dh = float(input_row["Domestic_Hot_Water_kW"])  # per end-user
    T_s = float(input_row["Supply_Temperature_degC"])
    T_r_sh = float(input_row["Return_Temperature_sh_degC"])
    T_r_dh = float(input_row["Return_Temperature_dhw_degC"])
    pump_head_bar = float(input_row["Pump_Head_Lift_bar"])
    end_user_dp_bar = float(input_row["End-User_Differential_Pressure_bar"])  # reserved but not used directly here

    # Mass flow per end-user at peak (kg/s)
    dh_sh = (h_water_kJ_per_kg(T_s) - h_water_kJ_per_kg(T_r_sh))  # kJ/kg
    dh_dh = (h_water_kJ_per_kg(T_s) - h_water_kJ_per_kg(T_r_dh))
    if dh_sh <= 0 or dh_dh <= 0:
        raise ValueError("Invalid temperature scheme leading to non-positive enthalpy difference.")
    m_dot_sh = (Q_sh) / (dh_sh / 1.0)  # kW / (kJ/kg) -> kg/s (since 1 kW = 1 kJ/s)
    m_dot_dh = (Q_dh) / (dh_dh / 1.0)

    # Pre-process network
    pipes_local = pipes.reset_index(drop=True).copy()
    source, leaves = find_source_and_leaves(pipes_local)
    pred_map = build_pred_map(pipes_local)

    # Critical route length
    lengths = pipes_local["length_m"].tolist()
    route_L = route_lengths_to_source(leaves, source, pred_map, lengths)
    L_critical = max(route_L)

    # Permissible pressure gradient (bar/m) along network
    max_grad_bar_per_m = (pump_head_bar - end_user_dp_bar) / L_critical

    # Consumers per node from service list
    # Sum ref_build per node
    node_consumers = (
        service.groupby("node_connc")["ref_build"].sum().astype(int).to_dict()
    )

    # Cumulative consumers carried by each edge (summed across routes)
    n_edges = len(pipes_local)
    CN_per_edge = cumulative_consumers_per_edge(leaves, source, pred_map, node_consumers, n_edges)

    # Mass flow per edge (kg/s)
    m_dot_edge: List[float] = []
    for CN in CN_per_edge:
        if CN <= 0:
            m_dot_edge.append(0.0)
            continue
        sf1 = sf_sh(CN)
        sf2 = sf_dhw(CN, substation_type)
        m_dot_edge.append(m_dot_sh * sf1 * CN + m_dot_dh * sf2 * CN)

    # Choose diameter per edge
    diam_col = []
    for (L_m, m_dot) in zip(pipes_local["length_m"].tolist(), m_dot_edge):
        sel_d = None
        for d_mm, k_mm in zip(catalogue["InnerDiameter_mm"], catalogue["Roughness_mm"]):
            dp_bar = pressure_loss_bar(L_m, float(d_mm), float(m_dot), float(k_mm), iters=friction_iters)
            grad = dp_bar / L_m if L_m > 0 else 0.0
            if grad <= max_grad_bar_per_m:
                sel_d = float(d_mm)
                break
        if sel_d is None:
            sel_d = float(catalogue["InnerDiameter_mm"].max())  # fallback: largest size
        diam_col.append(sel_d)

    pipe_results = pipes_local.copy()
    pipe_results["Selected_D_mm"] = diam_col

    # Dimension each service pipe independently (separate dimensioning as requested)
    s_diam: List[float] = []
    for _, row in service.iterrows():
        CN = int(row["ref_build"])  # no. of reference households at that service point
        L_m = float(row["length_m"])  # length of that service pipe
        if CN <= 0:
            s_diam.append(float(catalogue["InnerDiameter_mm"].min()))
            continue
        m_dot = m_dot_sh * sf_sh(CN) * CN + m_dot_dh * sf_dhw(CN, substation_type) * CN
        sel_d = None
        for d_mm, k_mm in zip(catalogue["InnerDiameter_mm"], catalogue["Roughness_mm"]):
            dp_bar = pressure_loss_bar(L_m, float(d_mm), float(m_dot), float(k_mm), iters=friction_iters)
            grad = dp_bar / L_m if L_m > 0 else 0.0
            if grad <= max_grad_bar_per_m:
                sel_d = float(d_mm)
                break
        if sel_d is None:
            sel_d = float(catalogue["InnerDiameter_mm"].max())
        s_diam.append(sel_d)

    service_results = service.copy()
    service_results["Selected_D_mm"] = s_diam

    meta = {
        "critical_length_m": float(L_critical),
        "max_pressure_gradient_bar_per_m": float(max_grad_bar_per_m),
        "source_node": int(source),
        "n_leaves": int(len(leaves)),
    }

    return pipe_results, service_results, meta


# -----------------------------
# Orchestration
# -----------------------------

def main(input_dir: Path = Path('.'), output_xlsx: Path = Path('dimensioning_results.xlsx'),
         substation_type: int = 1, friction_iters: int = 2) -> None:
    """Run the full workflow and write results to Excel with separate sheets."""
    input_data = read_input_data(input_dir / 'input-data.csv')
    pipes = load_pipe_segments(input_dir / 'pipe_segments.csv')
    nodes = load_nodes(input_dir / 'nodes.csv')
    service = load_service_pipes(input_dir / 'service_pipes.csv')
    catalogue = load_pipe_catalogue(optional_path=(input_dir / 'pipe_catalogue.csv'))

    pipe_res, serv_res, meta = dimension_network(
        input_data, pipes, nodes, service, catalogue,
        substation_type=substation_type, friction_iters=friction_iters)

    # Write to Excel (separate sheets)
    with pd.ExcelWriter(output_xlsx, engine='openpyxl', mode='w') as writer:
        # Input data as single-row DataFrame for clarity
        input_data.to_frame().T.to_excel(writer, sheet_name='InputData', index=False)
        pipe_res.to_excel(writer, sheet_name='PipeList', index=False)
        nodes.to_excel(writer, sheet_name='NodeList', index=False)
        serv_res.to_excel(writer, sheet_name='ServicePipes', index=False)
        # Meta-info sheet
        pd.DataFrame([meta]).to_excel(writer, sheet_name='Meta', index=False)

    print(f"Saved: {output_xlsx}")


if __name__ == '__main__':
    main()

