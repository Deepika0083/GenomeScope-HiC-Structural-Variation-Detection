"""
GenomeScope — Python SV Detection Script
=========================================
Project: Computational Framework for Detecting Structural Variations
         from 3D Genome (Hi-C) Data

This script does the same thing the browser does, but in pure Python.
You can run it on real Hi-C contact matrices (CSV format).

Usage:
    python sv_detector.py                         # runs on built-in demo data
    python sv_detector.py --matrix my_matrix.csv  # runs on your file
    python sv_detector.py --matrix my_matrix.csv --output results.csv

Requirements:
    pip install numpy scipy scikit-learn matplotlib pandas
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # works without a display (VS Code, servers)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from sklearn.ensemble import IsolationForest
import argparse
import os
import sys


# ─── STEP 1: Generate or load a contact matrix ─────────────────────────────

def make_demo_matrix(n=60, sv_type='del'):
    """
    Creates a synthetic Hi-C contact matrix with one planted SV.
    In a real project you'd load this from a .cool or .csv file.

    The matrix encodes how often each pair of genomic bins (10 kb each)
    was found in physical contact during the Hi-C experiment.
    """
    np.random.seed(42)
    matrix = np.zeros((n, n))

    # Normal distance-decay: nearby regions contact each other more
    for i in range(n):
        for j in range(i, n):
            dist = j - i
            val  = max(0, 1 - dist * 0.06) * np.random.lognormal(0, 0.3)
            matrix[i, j] = val
            matrix[j, i] = val

    # Plant an SV signal
    if sv_type == 'del':
        # Deletion: bright block off the diagonal (breakpoint fusion)
        matrix[10:18, 40:48] += 5.0
        matrix[40:48, 10:18] += 5.0

    elif sv_type == 'dup':
        # Duplication: parallel stripe
        matrix[8:15, 25:32] += 4.0
        matrix[25:32, 8:15] += 4.0

    elif sv_type == 'inv':
        # Inversion: anti-diagonal signal
        for k in range(n):
            mirror = n - 1 - k
            if 15 <= k <= 45 and 15 <= mirror <= 45:
                matrix[k, mirror] += 4.5
                matrix[mirror, k] += 4.5

    elif sv_type == 'tra':
        # Translocation: corner block (inter-chromosomal)
        matrix[5:12, 48:55] += 6.0
        matrix[48:55, 5:12] += 6.0

    return matrix


def load_matrix(path):
    """Load a contact matrix from a CSV or TSV file."""
    sep = '\t' if path.endswith('.tsv') else ','
    df  = pd.read_csv(path, header=None, sep=sep)
    return df.values.astype(float)


# ─── STEP 2: Normalize (observed / expected) ───────────────────────────────

def oe_normalize(matrix):
    """
    Observed / Expected normalization.

    For each genomic distance d, we compute the average contact frequency
    at that distance (the 'expected'), then divide every pixel by it.

    After this, the diagonal gradient disappears and real anomalies
    (SV breakpoints) stand out as bright patches.
    """
    n    = matrix.shape[0]
    oe   = np.zeros_like(matrix, dtype=float)

    for d in range(n):
        diag = np.diag(matrix, d)
        mean = diag.mean()
        if mean > 0:
            normed = diag / mean
            for idx, val in enumerate(normed):
                oe[idx, idx + d] = val
                oe[idx + d, idx] = val

    return oe


# ─── STEP 3: Anomaly detection ─────────────────────────────────────────────

def extract_window_features(oe_matrix, window=5):
    """
    Slide a window over the matrix and extract features for each pixel:
      - mean contact in the surrounding window
      - max contact
      - standard deviation
      - the pixel's own O/E value
      - z-score relative to the diagonal band
    """
    n        = oe_matrix.shape[0]
    features = []
    coords   = []
    W        = window

    for i in range(W, n - W):
        for j in range(i + 5, n - W):   # skip near-diagonal (normal TADs)
            patch = oe_matrix[i - W:i + W, j - W:j + W]
            # z-score: how many standard deviations above the mean at this distance?
            dist_band = np.diag(oe_matrix, j - i)
            z = (oe_matrix[i, j] - dist_band.mean()) / (dist_band.std() + 1e-9)

            features.append([
                patch.mean(),
                patch.max(),
                patch.std(),
                oe_matrix[i, j],
                z
            ])
            coords.append((i, j))

    return np.array(features), coords


def run_isolation_forest(features, contamination=0.03):
    """
    Isolation Forest: an unsupervised anomaly detection model.

    It works by randomly partitioning the data. Points that are easy to
    isolate (need fewer partitions) are considered anomalies — those are
    our SV candidates.

    contamination = expected fraction of anomalies (3% here).
    """
    clf    = IsolationForest(n_estimators=200, contamination=contamination, random_state=42)
    clf.fit(features)
    scores = clf.decision_function(features)   # more negative = more anomalous
    labels = clf.predict(features)             # -1 = anomaly, +1 = normal
    return scores, labels


def zscore_scan(oe_matrix, threshold=3.5):
    """
    Simple z-score scan across the O/E matrix.
    Pixels more than `threshold` standard deviations above the mean
    at their genomic distance are flagged.
    """
    n       = oe_matrix.shape[0]
    anomaly = np.zeros_like(oe_matrix, dtype=bool)

    for d in range(5, n):   # skip near-diagonal
        band = np.diag(oe_matrix, d)
        if band.std() > 0:
            z = (band - band.mean()) / band.std()
            for idx, zval in enumerate(z):
                if zval > threshold:
                    anomaly[idx, idx + d] = True
                    anomaly[idx + d, idx] = True

    return anomaly


# ─── STEP 4: Classify SV type from block position/shape ───────────────────

def classify_sv_type(i, j, n):
    """
    Simple rule-based classifier based on block position in the matrix.

    In a real implementation this would be a trained Random Forest or CNN.
    For the mini-project, these geometric rules work well enough:
      - Far off-diagonal (same chromosome)  → deletion or duplication
      - Anti-diagonal (mirror position)     → inversion
      - Very far corner blocks              → translocation
    """
    dist      = abs(j - i)
    rel_dist  = dist / n      # relative to matrix size

    if rel_dist > 0.6:
        return 'TRA', 'Translocation — possible inter-chromosomal fusion'
    elif rel_dist > 0.4:
        return 'DEL', 'Deletion — breakpoint fusion between distant regions'
    elif rel_dist > 0.25:
        return 'DUP', 'Duplication — elevated contact block near diagonal'
    else:
        return 'INV', 'Inversion — mirrored contact pattern'


# ─── STEP 5: Merge nearby calls into SV events ─────────────────────────────

def cluster_anomalies(coords, scores, n, bin_kb=10):
    """
    Merge pixels within 5 bins of each other into a single SV call.
    Returns a list of SV dicts.
    """
    if not coords:
        return []

    # Sort by anomaly score (most anomalous first)
    sorted_hits = sorted(zip(scores, coords), key=lambda x: x[0])[:50]

    svs    = []
    used   = set()
    sv_id  = 1

    for score, (ci, cj) in sorted_hits:
        if (ci, cj) in used:
            continue

        # Collect nearby pixels
        cluster_i, cluster_j = [ci], [cj]
        for score2, (ci2, cj2) in sorted_hits:
            if abs(ci2 - ci) <= 5 and abs(cj2 - cj) <= 5:
                cluster_i.append(ci2)
                cluster_j.append(cj2)
                used.add((ci2, cj2))
        used.add((ci, cj))

        # Compute SV coordinates
        i_start = min(cluster_i)
        i_end   = max(cluster_i)
        j_start = min(cluster_j)
        j_end   = max(cluster_j)
        size_kb = (j_end - i_start) * bin_kb

        if size_kb < 50:   # skip very small calls
            continue

        sv_type, description = classify_sv_type((i_start + i_end) // 2,
                                                 (j_start + j_end) // 2, n)
        confidence = min(0.99, max(0.40, abs(score) * 3))

        svs.append({
            'SV_ID':      f'SV{sv_id:03d}',
            'Type':       sv_type,
            'Bin_Start':  i_start,
            'Bin_End':    j_end,
            'Start_kb':   i_start * bin_kb,
            'End_kb':     j_end   * bin_kb,
            'Size_kb':    size_kb,
            'Confidence': round(confidence, 2),
            'Note':       description,
        })
        sv_id += 1
        if sv_id > 10:
            break

    return svs


# ─── STEP 6: Visualize ─────────────────────────────────────────────────────

def plot_results(raw_matrix, oe_matrix, svs, output_path='genomescope_output.png'):
    """
    Saves a 3-panel figure:
      Left  – raw contact matrix
      Middle – O/E normalized matrix
      Right  – O/E with SV breakpoints marked
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor('#07090f')

    def draw_heatmap(ax, data, title, cmap='hot'):
        ax.set_facecolor('#07090f')
        im = ax.imshow(data, cmap=cmap, aspect='auto',
                       vmin=np.percentile(data, 5), vmax=np.percentile(data, 98))
        ax.set_title(title, color='#e8eaf6', fontsize=11, pad=8)
        ax.tick_params(colors='#4a5568')
        for spine in ax.spines.values():
            spine.set_edgecolor('#1a2035')
        return im

    draw_heatmap(axes[0], raw_matrix,  'Raw contact matrix')
    draw_heatmap(axes[1], oe_matrix,   'After O/E normalization')

    sv_overlay = oe_matrix.copy()
    ax3 = axes[2]
    draw_heatmap(ax3, sv_overlay, 'Detected SVs')

    # SV colour map
    type_colors = {'DEL': '#ef4444', 'DUP': '#4f9cf9', 'INV': '#f59e0b', 'TRA': '#22c55e'}

    for sv in svs:
        color = type_colors.get(sv['Type'], 'white')
        rect  = mpatches.Rectangle(
            (sv['Bin_Start'] - 0.5, sv['Bin_Start'] - 0.5),
            sv['Bin_End'] - sv['Bin_Start'],
            sv['Bin_End'] - sv['Bin_Start'],
            linewidth=2, edgecolor=color, facecolor='none', linestyle='--'
        )
        ax3.add_patch(rect)
        ax3.text(sv['Bin_Start'], sv['Bin_Start'] - 1,
                 f"{sv['SV_ID']} {sv['Type']}",
                 color=color, fontsize=7, fontfamily='monospace')

    # Legend
    patches = [mpatches.Patch(color=c, label=t) for t, c in type_colors.items()]
    axes[2].legend(handles=patches, loc='lower right',
                   facecolor='#0d1120', edgecolor='#1a2035',
                   labelcolor='#e8eaf6', fontsize=8)

    plt.tight_layout(pad=2)
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#07090f')
    print(f"\n  Plot saved → {output_path}")


# ─── STEP 7: Print and save report ─────────────────────────────────────────

def print_report(svs):
    print("\n" + "═" * 70)
    print("  GENOMESCOPE — SV Detection Report")
    print("═" * 70)
    if not svs:
        print("  No SVs detected above threshold.")
        return

    header = f"  {'SV_ID':<8} {'Type':<5} {'Start(kb)':<12} {'End(kb)':<10} {'Size(kb)':<10} {'Conf':>6}  Note"
    print(header)
    print("  " + "-" * 68)
    for sv in svs:
        print(f"  {sv['SV_ID']:<8} {sv['Type']:<5} {sv['Start_kb']:<12} {sv['End_kb']:<10} {sv['Size_kb']:<10} {sv['Confidence']:>6.0%}  {sv['Note'][:40]}")
    print("═" * 70)


def save_csv(svs, path):
    if not svs:
        return
    df = pd.DataFrame(svs)[['SV_ID', 'Type', 'Start_kb', 'End_kb', 'Size_kb', 'Confidence', 'Note']]
    df.to_csv(path, index=False)
    print(f"  CSV saved → {path}")


# ─── MAIN ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='GenomeScope SV Detector')
    parser.add_argument('--matrix',       default=None,               help='Path to contact matrix CSV')
    parser.add_argument('--sv_type',      default='del',              help='Demo SV type: del/dup/inv/tra (ignored if --matrix given)')
    parser.add_argument('--output',       default='sv_report.csv',    help='Output CSV path')
    parser.add_argument('--plot',         default='genomescope_output.png', help='Output plot path')
    parser.add_argument('--sensitivity',  default=3.5, type=float,    help='Z-score threshold (default 3.5)')
    args = parser.parse_args()

    print("\n  ┌─────────────────────────────────────────┐")
    print("  │  GenomeScope — SV Detection Pipeline    │")
    print("  └─────────────────────────────────────────┘\n")

    # 1. Load or generate matrix
    if args.matrix:
        print(f"  [1/6] Loading matrix from {args.matrix}...")
        matrix = load_matrix(args.matrix)
    else:
        print(f"  [1/6] Generating demo matrix (SV type: {args.sv_type})...")
        matrix = make_demo_matrix(n=60, sv_type=args.sv_type)
    n = matrix.shape[0]
    print(f"        Matrix size: {n}×{n} bins")

    # 2. Normalize
    print("  [2/6] Applying O/E normalization...")
    oe = oe_normalize(matrix)

    # 3. Feature extraction
    print("  [3/6] Extracting window features...")
    features, coords = extract_window_features(oe, window=4)
    print(f"        {len(features)} candidate windows")

    # 4. Anomaly detection
    print("  [4/6] Running Isolation Forest...")
    if_scores, if_labels = run_isolation_forest(features, contamination=0.04)

    print("  [4/6] Running z-score scan...")
    z_anomaly = zscore_scan(oe, threshold=args.sensitivity)

    # Combine: use IF anomaly coords
    anomaly_coords = [coords[k] for k in range(len(if_labels)) if if_labels[k] == -1]
    anomaly_scores = [if_scores[k] for k in range(len(if_labels)) if if_labels[k] == -1]

    # 5. Cluster into SVs
    print("  [5/6] Clustering anomalies into SV calls...")
    svs = cluster_anomalies(anomaly_coords, anomaly_scores, n)

    # 6. Output
    print("  [6/6] Generating output...\n")
    print_report(svs)
    save_csv(svs, args.output)
    plot_results(matrix, oe, svs, args.plot)

    print(f"\n  Done. {len(svs)} SV(s) detected.\n")


if __name__ == '__main__':
    main()
