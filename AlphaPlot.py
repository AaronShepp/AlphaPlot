import random # randomness (shuffling, swaps, probabilities)
import os # file path management
import math # math functions (exp, sqrt)
import time # timing
import pandas as pd # tabular data export
import matplotlib.pyplot as plt # plotting

# ==========================================================
# ========================= INPUTS =========================
# ==========================================================

# ---------- Geometry ----------
INNER_DIAMETER = 9 #must be odd, minimum 3
INNER_RADIUS   = (INNER_DIAMETER - 1) // 2
EDGE_WIDTH     = 3 #optional

# ---------- Species list ----------
LISTS = {
"SPECIES_LIST_1": ["Sorbus torminalis", "Taxus baccata", "Acer campestre", "Tilia platyphyllos"],
"SPECIES_LIST_2": ["Sorbus torminalis", "Taxus baccata", "Pinus nigra", "Prunus avium"],
"SPECIES_LIST_3": ["Acer campestre", "Tilia platyphyllos", "Pinus nigra", "Prunus avium"],
"SPECIES_LIST_4": ["Sorbus torminalis", "Taxus baccata", "Acer campestre", "Tilia platyphyllos", "Pinus nigra", "Prunus avium"]
}

CHOOSE_LIST = "SPECIES_LIST_1"
LIST = LISTS[CHOOSE_LIST]

# ---------- Plot colors ----------
SPECIES_COLORS = {
    "Sorbus torminalis":  "#ff0000",
    "Taxus baccata":   "#020202",
    "Acer campestre":    "#ffff00",
    "Tilia platyphyllos":   "#16f138",
    "Pinus nigra":   "#0000ff",
    "Prunus avium":  "#e816f1"
}

# ---------- Criteria ----------
MIN_NEIGHBOR_CASES = {0: 1, #for n = 0, mink >= _
                      1: 1, #for n = 1, mink >= _
                      2: 1, 
                      3: 0, 
                      4: 0,
                      5: 0,
                      6: 0}

NEIGHBOR_TOLERANCE = {0: 4, #for n = 0, rangek <= _
                      1: 4, 
                      2: 4, 
                      3: 2,
                      4: 1,
                      5: 1,
                      6: 0}

# ---------- Scoring ----------
PERFECT_BONUS                 = 500
POINTS_PER_MIN_CONSTRAINT     = 20 # per FocalSpecies-NeighbborSpecies-minK_Success. Does not incentivize over-optimization
POINTS_PER_BALANCE_CONSTRAINT = 10 # per FocalSpecies-NeighbborSpecies-rangeK_Success. Multiplied by how much the criteria is surpassed +1

MIN_VIOLATION_PENALTY = 200 # per FocalSpecies-NeighbborSpecies-minK_Violation. Multiplied by how much the criteria is violated
TOL_VIOLATION_PENALTY = 20 # per FocalSpecies-NeighbborSpecies-rangeK_Violation. Multiplied by how much the criteria is violated.

VARIANCE_PENALTY = 2.0 # across all FocalSpecies-NeighbborSpecies case counts
RANGE_PENALTY    = 5.0 # penalizes high range case counts for a given n
GEN_IMPROVE_WEIGHT   = 2.0 # rewards high minima of case counts (raising worst-off k for a given n), prioritzes n with low minima
MAX_EVALUATED_NEIGHBORS = 3 # Used for ecological weight ie. "raise the weakest case". Not used for MIN and TOLERANCE criteria

# ---------- Simulated Annealing Calibration ----------
# Temperature(T): likeliness to accept worse swaps, Probobility ​= e^(Δ/T), cooling is scheduled (geometric) 
SA_T0    =  1200.0 # initial temperature, recommended: equivalent to average score difference between random swaps
SA_T_END = 0.1 #minimum temperature

# ---------- Search Effort ----------
RANDOM_ATTEMPTS = 1000
SA_SWAPS        = 5000
SEED_RANGE = range(2)

# ---------- Export ----------
TEST = True # If True, adds 'Test' prefix to output filenames.

OUTPUT_NAME = (
    f"{'Test' if TEST else ''}"
    f"{CHOOSE_LIST}_"
    f"c{INNER_DIAMETER}"
    f"e{EDGE_WIDTH}"
    f"a{RANDOM_ATTEMPTS}"
    f"sw{SA_SWAPS}"
    f"s{len(SEED_RANGE)}"
)
#Unique filenames created on conflict

EXPORT_DIR = r"C:\MyExperiment\Method\PlantingLayout\AlphaPlot\Exports\c9e3\Tests" #Where do you want the exports to be saved? 

#Remember to close pop-up window of exported plot or you will not be able to run it again.
#If no plot that fits criteria is found, no exports will be made.

# ==========================================================
# ===================== CONSTANTS ==========================
# ==========================================================

NEIGHBOR_OFFSETS = [ # Define axial neighbor offsets
    (+1, 0), (+1, -1), (0, -1),
    (-1, 0), (-1, +1), (0, +1)
]

# ==========================================================
# ====================== GEOMETRY ==========================
# ==========================================================

def generate_hexagon(radius): # Creates all axial coordinates inside a hexagon of given radius.
    coords = set() # Stores coordinates without duplicates.
    for q in range(-radius, radius + 1): # Iterate over axial q-coordinates.
        r1 = max(-radius, -q - radius)
        r2 = min(radius, -q + radius)
        for r in range(r1, r2 + 1):
            coords.add((q, r))
    return coords


def split_inner_and_edge(inner_radius, edge_width):
    inner = generate_hexagon(inner_radius)
    total = generate_hexagon(inner_radius + edge_width)
    edge  = total - inner
    return inner, edge, total


def neighbors(coord):
    q, r = coord
    return [(q + dq, r + dr) for dq, dr in NEIGHBOR_OFFSETS]

def axial_distance(q, r):
    return max(abs(q), abs(r), abs(-q - r))


def axial_spiral_order(coords):
    """
    Return coords ordered in an outward hex spiral,
    starting at (0, 0), clockwise, ring by ring.
    """
    coords = set(coords)
    ordered = []

    if (0, 0) in coords:
        ordered.append((0, 0))

    # Directions for walking a hex ring (clockwise)
    directions = [
        (0, -1),
        (-1, 0),
        (-1, 1),
        (0, 1),
        (1, 0),
        (1, -1),
    ]

    max_radius = max(axial_distance(q, r) for q, r in coords)

    for radius in range(1, max_radius + 1):

        # Start at "east-north-east" corner of the ring
        q, r = radius, 0

        for dq, dr in directions:
            for _ in range(radius):
                if (q, r) in coords:
                    ordered.append((q, r))
                q += dq
                r += dr

    return ordered


# ==========================================================
# ===================== ASSIGNMENT =========================
# ==========================================================

def balanced_assignment(coords, species_list):
    coords = list(coords)
    reps, rem = divmod(len(coords), len(species_list))
    pool = [sp for sp in species_list for _ in range(reps)]
    pool += species_list[:rem]
    random.shuffle(pool)
    return dict(zip(coords, pool))

# ==========================================================
# ================= NEIGHBOR COUNTING ======================
# ==========================================================

def compute_adjacency(inner_coords, species_map, species_list):
    counts = {
        sp: {other: {n_neighbors: 0 for n_neighbors in range(7)} for other in species_list}
        for sp in species_list
    }

    for coord in inner_coords:
        focal = species_map[coord]
        neigh = [species_map[n] for n in neighbors(coord) if n in species_map]
        for other in species_list:
            counts[focal][other][neigh.count(other)] += 1

    return counts

# ==========================================================
# ==================== HARD CONSTRAINT =====================
# ==========================================================

def satisfies_constraints(counts, species_list):
    for sp in species_list:
        for other in species_list:
            for n_neighbors, m in MIN_NEIGHBOR_CASES.items():
                if counts[sp][other][n_neighbors] < m:
                    return False

    for other in species_list:
        for n_neighbors, tol in NEIGHBOR_TOLERANCE.items():
            vals = [counts[sp][other][n_neighbors] for sp in species_list]
            if max(vals) - min(vals) > tol:
                return False

    return True

# ==========================================================
# ======================= SCORING ==========================
# ==========================================================

def score_layout(counts, species_list):
    score = 0
    all_satisfied = True

    for sp in species_list:
        for other in species_list:
            for n_neighbors, m in MIN_NEIGHBOR_CASES.items():
                delta = counts[sp][other][n_neighbors] - m
                if delta >= 0:
                    score += POINTS_PER_MIN_CONSTRAINT
                else:
                    score -= MIN_VIOLATION_PENALTY * abs(delta)
                    all_satisfied = False

    for other in species_list:
        for n_neighbors, tol in NEIGHBOR_TOLERANCE.items():
            vals = [counts[sp][other][n_neighbors] for sp in species_list]
            r = max(vals) - min(vals)

            if r <= tol:
                score += POINTS_PER_BALANCE_CONSTRAINT * (tol - r + 1)
            else:
                score -= TOL_VIOLATION_PENALTY * (r - tol)
                all_satisfied = False

            mean = sum(vals) / len(vals)
            variance = sum((v - mean) ** 2 for v in vals)
            score -= VARIANCE_PENALTY * variance
            score -= RANGE_PENALTY * r

    for other in species_list:
        mins = {n_neighbors: min(counts[sp][other][n_neighbors] for sp in species_list)
                for n_neighbors in range(MAX_EVALUATED_NEIGHBORS + 1)}
        max_min = max(mins.values())
        for n_neighbors, m in mins.items():
            score += GEN_IMPROVE_WEIGHT * (10 if m < max_min else 2) * m

    if all_satisfied:
        score += PERFECT_BONUS

    return score

# ==========================================================
# ================= SIMULATED ANNEALING ====================
# ==========================================================

def random_swap(species_map, inner_coords):
    a, b = random.sample(inner_coords, 2)
    new = species_map.copy()
    new[a], new[b] = new[b], new[a]
    return new


def simulated_annealing(start_map, inner_coords, species_list):
    current = start_map
    current_score = score_layout(
        compute_adjacency(inner_coords, current, species_list),
        species_list
    )

    best = current
    best_score = current_score
    inner_coords = list(inner_coords)

    for step in range(SA_SWAPS):
        T = SA_T0 * (SA_T_END / SA_T0) ** (step / SA_SWAPS)
        candidate = random_swap(current, inner_coords)
        counts = compute_adjacency(inner_coords, candidate, species_list)
        cand_score = score_layout(counts, species_list)

        if cand_score > current_score or random.random() < math.exp((cand_score - current_score) / T):
            current, current_score = candidate, cand_score

        if cand_score > best_score:
            best, best_score = candidate, cand_score

    return best, best_score

# ==========================================================
# ==================== SEARCH PIPELINE =====================
# ==========================================================

def find_best_layout(inner_coords, edge_coords, species_list):
    best_overall = (-float("inf"), None)
    best_valid   = (-float("inf"), None)

    for _ in range(RANDOM_ATTEMPTS):
        layout = {
            **balanced_assignment(inner_coords, species_list),
            **balanced_assignment(edge_coords, species_list)
        }
        counts = compute_adjacency(inner_coords, layout, species_list)
        score = score_layout(counts, species_list)
        valid = satisfies_constraints(counts, species_list)

        if score > best_overall[0]:
            best_overall = (score, layout)

        if valid and score > best_valid[0]:
            best_valid = (score, layout)

    start_layout = best_valid[1] if best_valid[1] is not None else best_overall[1]

    improved, improved_score = simulated_annealing(
        start_layout, inner_coords, species_list
    )

    final_counts = compute_adjacency(inner_coords, improved, species_list)
    return improved, final_counts, improved_score

# ==========================================================
# ================== NEIGHBOR DIVERSITY ====================
# ==========================================================
def shannon_index(counts):
    total = sum(counts.values())
    if total == 0:
        return 0
    H = -sum((c/total) * math.log(c/total) for c in counts.values() if c > 0)
    return H

def pielou_evenness(counts):
    richness = sum(1 for c in counts.values() if c > 0)
    if richness <= 1:
        return 0
    H = shannon_index(counts)
    return H / math.log(richness)


# ==========================================================
# ======================= DRAW PLOT ===========================
# ==========================================================

def avoid_overwrite(path):
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    i = 1
    while True:
        new = f"{base}_{i}{ext}"
        if not os.path.exists(new):
            return new
        i += 1


def axial_to_pixel(q, r):
    return (
        math.sqrt(3) * q + math.sqrt(3) / 2 * r,
        1.5 * r
    )

def hex_boundary_pixels(radius):
    """
    Return pixel-space vertices of a hexagon boundary at a given axial radius.
    Radius may be fractional (e.g. N + 0.5).
    """
    corners = [
        ( radius,  0),
        ( radius, -radius),
        ( 0,      -radius),
        (-radius,  0),
        (-radius,  radius),
        ( 0,       radius),
    ]
    return [axial_to_pixel(q, r) for q, r in corners]


from matplotlib.patches import Polygon

def plot_layout(coords, species_map):
    fig, ax = plt.subplots(figsize=(12, 12))

    # ---- Correct hex boundaries (between cells) ----
    inner_boundary = INNER_RADIUS + 0.5
    outer_boundary = INNER_RADIUS + EDGE_WIDTH + 0.5

    # Outer grey ring with black edge only on outer boundary
    outer_hex = hex_boundary_pixels(outer_boundary)
    inner_hex = hex_boundary_pixels(inner_boundary)

    # Grey fill
    ax.add_patch(
        Polygon(
            outer_hex,
            closed=True,
            facecolor="#ababab",
            edgecolor="none",  # no edge here
            zorder=0
        )
    )

    # White inner core
    ax.add_patch(
        Polygon(
            inner_hex,
            closed=True,
            facecolor="white",
            edgecolor="none",
            zorder=1
        )
    )

    # Black line along **outer boundary**
    ax.add_patch(
        Polygon(
            outer_hex,
            closed=True,
            fill=False,       # no fill, just the edge
            edgecolor="black",
            linewidth=1,
            zorder=2
        )
    )

    
    # ---- Scatter trees ----
    marker_size = -40 * (INNER_DIAMETER + 2 * EDGE_WIDTH) + 880
    marker_size = max(marker_size, 100)  # optional safety floor

    xs, ys, colors = [], [], []
    for q, r in coords:
        x, y = axial_to_pixel(q, r)
        xs.append(x)
        ys.append(y)
        colors.append(SPECIES_COLORS[species_map[(q, r)]])

    ax.scatter(
        xs, ys,
        s=marker_size,
        c=colors,
        edgecolors="k",
        linewidths=0.5,
        zorder=3
    )

    ax.set_aspect("equal")
    ax.axis("off")


# ==========================================================
# ========================= MAIN ===========================
# ==========================================================

if __name__ == "__main__":

    start = time.time()
    os.makedirs(EXPORT_DIR, exist_ok=True)

    inner, edge, all_coords = split_inner_and_edge(INNER_RADIUS, EDGE_WIDTH)

    global_best_valid   = (-float("inf"), None, None)
    global_best_overall = (-float("inf"), None, None)

    for seed in SEED_RANGE:
        random.seed(seed)
        layout, counts, score = find_best_layout(inner, edge, LIST)
        valid = satisfies_constraints(counts, LIST)

        print(f"Seed {seed}: score={score:.2f}, valid={valid}")

        # Update best valid layout
        if valid and score > global_best_valid[0]:
            global_best_valid = (score, layout, counts)

        # Always update best overall layout
        if score > global_best_overall[0]:
            global_best_overall = (score, layout, counts)

    if global_best_valid[1] is not None:
        best_score, final_map, final_counts = global_best_valid
        print("Exporting BEST VALID layout")
    else:
        best_score, final_map, final_counts = global_best_overall
        print("Exporting BEST OVERALL layout (criteria not fully met)")

    if final_map is None and global_best_overall[1] is not None:
        print("No layout fully satisfies constraints. Exporting best overall layout.")
        best_score, final_map, final_counts = global_best_overall
    elif final_map is None:
        raise RuntimeError("No layouts could be generated at all.")


    print("\n=== FINAL EXPORTED LAYOUT ===")
    if satisfies_constraints(final_counts, LIST):
        print("Exporting BEST VALID layout (all constraints met)")
    else:
        print("Exporting BEST OVERALL layout (constraints not fully met)")
    print(f"Score = {best_score:.2f}")
    print(f"Completed in {time.time() - start:.2f} seconds")
    
    img_path   = avoid_overwrite(os.path.join(EXPORT_DIR, OUTPUT_NAME + ".svg"))
    plot_layout(all_coords, final_map)
    plt.savefig(img_path, format="svg", transparent=True)    

    # ================== Criteria Statistics ==============
     
    rows = []
    for sp in LIST:
        for other in LIST:
            row = {"Focal": sp, "Neighbor": other}
            for n_neighbors in range(7):
                row[f"{n_neighbors}_neighbors"] = final_counts[sp][other][n_neighbors]
            rows.append(row)

    df = pd.DataFrame(rows)

    excel_path = avoid_overwrite(os.path.join(EXPORT_DIR, OUTPUT_NAME + "_CritStats.xlsx"))
    import numpy as np

    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:

        df.to_excel(writer, index=False, sheet_name="NeighborCounts")

        workbook  = writer.book
        worksheet = writer.sheets["NeighborCounts"]

        # ----- Formatting styles -----
        bold = workbook.add_format({"bold": True})
        yellow      = workbook.add_format({"bg_color": "#fff000", "bold": True})
        green = workbook.add_format({"bg_color": "#92d050", "bold": True})
        blue  = workbook.add_format({"bg_color": "#00b0f0", "bold": True})

        start_row = len(df) + 1

        # ----- Write criteria rows -----
        worksheet.write(start_row, 0, "MIN_CRITERIA", bold)
        worksheet.write(start_row + 1, 0, "RANGE_CRITERIA", bold)
        worksheet.write(start_row + 2, 0, "MIN", bold)
        worksheet.write(start_row + 3, 0, "RANGE", bold)

        for n in range(7):

            col = 2 + n
            column_values = df.iloc[:, col].values

            min_crit  = MIN_NEIGHBOR_CASES[n]
            range_crit = NEIGHBOR_TOLERANCE[n]

            col_min   = column_values.min()
            col_range = column_values.max() - column_values.min()

            # ---- Write criteria ----
            worksheet.write(start_row,     col, min_crit, bold)
            worksheet.write(start_row + 1, col, range_crit, bold)

            # ---- Choose MIN formatting ----
            if col_min < min_crit:
                fmt_min = yellow
            elif col_min == min_crit:
                fmt_min = green
            else:
                fmt_min = blue

            worksheet.write(start_row + 2, col, col_min, fmt_min)

            # ---- Choose RANGE formatting ----
            if col_range > range_crit:
                fmt_range = yellow
            elif col_range == range_crit:
                fmt_range = green
            else:
                fmt_range = blue

            worksheet.write(start_row + 3, col, col_range, fmt_range)
    
    # ================== Cell Statistics ==============
    cell_rows = []

    ordered_coords = axial_spiral_order(final_map.keys())

    for coord in ordered_coords:
        sp = final_map[coord]

        neigh = [final_map[n] for n in neighbors(coord) if n in final_map]
        counts = {s: neigh.count(s) for s in LIST}

        richness = sum(1 for v in counts.values() if v > 0)
        
        def dominant_neighbor(neighbor_counts):
            """
            neighbor_counts: dict {species_name: count}

            Returns:
                str: dominant species, or comma-separated ties
            """
            if not neighbor_counts:
                return ""

            max_val = max(neighbor_counts.values())

            tied = [
                species for species, count in neighbor_counts.items()
                if count == max_val
            ]

            return ", ".join(tied)

        dominant = dominant_neighbor(counts)


        cell_rows.append({
            "q": coord[0],
            "r": coord[1],
            "species": sp,
            "zone": "core" if coord in inner else "edge",
            "hex_radius": axial_distance(coord[0], coord[1]),
            "neighbor_richness": richness,
            "shannon_index": shannon_index(counts),
            "pielou_evenness": pielou_evenness(counts),
            "dominant_neighbor_species": dominant
        })

    cell_df = pd.DataFrame(cell_rows)
    cell_path = avoid_overwrite(os.path.join(EXPORT_DIR, OUTPUT_NAME + "_Layout.xlsx"))
    cell_df.to_excel(cell_path, index=False)

# ================== Parameter Record ==============

    txt_path = avoid_overwrite(os.path.join(EXPORT_DIR, OUTPUT_NAME + "_Params.txt"))

    with open(txt_path, "w") as f:
        f.write("=== Geometry ===\n")
        f.write(f"INNER_DIAMETER = {INNER_DIAMETER}\n")
        f.write(f"EDGE_WIDTH     = {EDGE_WIDTH}\n\n")

        f.write("=== Species list ===\n")
        f.write(f"LIST = {LIST}\n\n")

        f.write("=== Criteria ===\n")
        f.write(f"MIN_NEIGHBOR_CASES = { {k: MIN_NEIGHBOR_CASES[k] for k in range(1, 7)} }\n")
        f.write(f"NEIGHBOR_TOLERANCE = { {k: NEIGHBOR_TOLERANCE[k] for k in range(1, 7)} }\n\n")

        f.write("=== Scoring ===\n")
        f.write(f"PERFECT_BONUS                 = {PERFECT_BONUS}\n")
        f.write(f"POINTS_PER_MIN_CONSTRAINT     = {POINTS_PER_MIN_CONSTRAINT}\n")
        f.write(f"POINTS_PER_BALANCE_CONSTRAINT = {POINTS_PER_BALANCE_CONSTRAINT}\n")
        f.write(f"MIN_VIOLATION_PENALTY         = {MIN_VIOLATION_PENALTY}\n")
        f.write(f"TOL_VIOLATION_PENALTY         = {TOL_VIOLATION_PENALTY}\n")
        f.write(f"VARIANCE_PENALTY              = {VARIANCE_PENALTY}\n")
        f.write(f"RANGE_PENALTY                 = {RANGE_PENALTY}\n")
        f.write(f"GEN_IMPROVE_WEIGHT            = {GEN_IMPROVE_WEIGHT}\n")
        f.write(f"MAX_EVALUATED_NEIGHBORS       = {MAX_EVALUATED_NEIGHBORS}\n\n")

        f.write("=== Simulated Annealing Calibration ===\n")
        f.write(f"SA_T0    = {SA_T0}\n")
        f.write(f"SA_T_END = {SA_T_END}\n\n")

        f.write("=== Search Effort ===\n")
        f.write(f"RANDOM_ATTEMPTS = {RANDOM_ATTEMPTS}\n")
        f.write(f"SA_SWAPS        = {SA_SWAPS}\n")
        f.write(f"SEED_RANGE      = {list(SEED_RANGE)}\n")

    
    print(f"Exports Path:  {img_path}")
  
    plt.show()

