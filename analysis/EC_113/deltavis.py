#!/usr/bin/env python3
"""DeltaVis: matplotlib dot plot for nucmer/MUMmer .delta alignment files.

Usage:
    python deltavis.py input.delta [--rotate OFFSET] [--flip] [-o output.png]

The --rotate flag specifies a circular coordinate rotation offset (in bp)
for the query (y-axis). This shifts all query coordinates by the given
amount, wrapping around the total query length.

The --flip flag reverses the query axis, as if the reverse complement
had been aligned. Tick labels still show original coordinates.
"""

import argparse
import sys
from collections import OrderedDict
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class Alignment:
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    errors: int
    sim_errors: int
    stop_codons: int


@dataclass
class AlignmentSection:
    ref_id: str
    query_id: str
    ref_len: int
    query_len: int
    alignments: list = field(default_factory=list)


@dataclass
class DeltaFile:
    reference_path: str
    query_path: str
    fmt: str
    sections: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_delta(path: str) -> DeltaFile:
    with open(path) as f:
        lines = f.read().splitlines()

    if len(lines) < 2:
        raise ValueError("Delta file too short")

    paths = lines[0].split()
    ref_path, query_path = paths[0], paths[-1]

    fmt = lines[1].strip()
    if fmt not in ("NUCMER", "PROMER"):
        raise ValueError(f"Expected NUCMER or PROMER, got '{fmt}'")

    delta = DeltaFile(ref_path, query_path, fmt)
    current_section = None
    in_distances = False

    for line in lines[2:]:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Flush previous section
            if current_section is not None:
                delta.sections.append(current_section)

            parts = line[1:].split()
            ref_id, query_id = parts[0], parts[1]
            ref_len, query_len = int(parts[2]), int(parts[3])
            current_section = AlignmentSection(ref_id, query_id, ref_len, query_len)
            in_distances = False

        elif in_distances:
            val = int(line)
            if val == 0:
                in_distances = False

        else:
            parts = line.split()
            if len(parts) == 7:
                ints = list(map(int, parts))
                current_section.alignments.append(Alignment(*ints))
                in_distances = True
            elif len(parts) == 1:
                val = int(parts[0])
                if val == 0:
                    in_distances = False
            else:
                raise ValueError(f"Unexpected line: {line}")

    if current_section is not None:
        delta.sections.append(current_section)

    return delta


# ---------------------------------------------------------------------------
# Layout computation
# ---------------------------------------------------------------------------

def compute_layout(delta: DeltaFile):
    """Return (refs, queries, ref_offsets, query_offsets, total_ref, total_query).

    refs/queries: OrderedDict name -> length
    *_offsets: dict name -> cumulative offset
    """
    refs = OrderedDict()
    queries = OrderedDict()

    for sec in delta.sections:
        if sec.ref_id not in refs:
            refs[sec.ref_id] = sec.ref_len
        if sec.query_id not in queries:
            queries[sec.query_id] = sec.query_len

    def cumulative(seqs):
        offsets = {}
        total = 0
        for name, length in seqs.items():
            offsets[name] = total
            total += length
        return offsets, total

    ref_offsets, total_ref = cumulative(refs)
    query_offsets, total_query = cumulative(queries)

    return refs, queries, ref_offsets, query_offsets, total_ref, total_query


# ---------------------------------------------------------------------------
# Auto-detection of flip and rotation
# ---------------------------------------------------------------------------

def auto_orient(delta: DeltaFile):
    """Detect optimal flip and rotation for the query (y-axis).

    Returns (rotate, flip) where rotate is in bp and flip is a bool.

    Strategy:
    - Flip: choose the orientation (forward vs reverse) with more total
      aligned bases.
    - Rotation: maximize a monotonicity score — the fraction of consecutive
      main-diagonal alignment pairs (sorted by reference position, weighted
      by alignment length) whose display-y values are non-decreasing.
    """
    refs, queries, ref_offsets, query_offsets, total_ref, total_query = compute_layout(delta)
    if total_ref == 0 or total_query == 0:
        return 0, False

    scale = total_query / total_ref

    # Gather per-alignment info
    fwd_bases = 0
    rev_bases = 0
    all_info = []  # (diag_offset, ref_mid, query_mid, aln_len)
    for sec in delta.sections:
        r_off = ref_offsets[sec.ref_id]
        q_off = query_offsets[sec.query_id]
        for aln in sec.alignments:
            q_start, q_end = aln.query_start, aln.query_end
            is_forward = q_start <= q_end
            aln_len = abs(aln.ref_end - aln.ref_start)
            if is_forward:
                fwd_bases += aln_len
            else:
                rev_bases += aln_len
            ref_mid = r_off + (aln.ref_start + aln.ref_end) / 2
            query_mid = q_off + (q_start + q_end) / 2
            all_info.append((ref_mid, query_mid, aln_len))

    flip = rev_bases > fwd_bases

    if not all_info:
        return 0, flip

    # Compute diagonal offsets to identify main-diagonal cluster
    diag_data = []
    for ref_mid, query_mid, aln_len in all_info:
        effective_q = (total_query - query_mid) if flip else query_mid
        off = (effective_q - ref_mid * scale) % total_query
        diag_data.append((off, ref_mid, query_mid, aln_len))

    n_bins = 36
    bin_width = total_query / n_bins
    hist = np.zeros(n_bins)
    for off, _, _, w in diag_data:
        hist[int(off / bin_width) % n_bins] += w
    peak_bin = int(np.argmax(hist))
    threshold = hist[peak_bin] * 0.1
    main_bins = {i for i in range(n_bins) if hist[i] >= threshold}

    # Keep only main-diagonal alignments, sorted by ref position
    main = sorted(
        [(rm, qm, w) for off, rm, qm, w in diag_data
         if int(off / bin_width) % n_bins in main_bins],
        key=lambda x: x[0])

    if len(main) < 2:
        return 0, flip

    # Pre-compute pair weights for scoring
    pair_weights = [min(main[i][2], main[i + 1][2])
                    for i in range(len(main) - 1)]
    total_weight = sum(pair_weights)
    if total_weight == 0:
        return 0, flip

    def monotonicity_score(R):
        displays = []
        for _, query_mid, _ in main:
            if flip:
                displays.append(total_query - (query_mid - R) % total_query)
            else:
                displays.append((query_mid - R) % total_query)
        score = 0.0
        for i in range(len(displays) - 1):
            if displays[i + 1] >= displays[i]:
                score += pair_weights[i]
        return score / total_weight

    # Coarse search
    n_coarse = 1000
    coarse_rots = np.linspace(0, total_query, n_coarse, endpoint=False)
    best_score = -1
    best_r = 0
    for r in coarse_rots:
        s = monotonicity_score(int(r))
        if s > best_score:
            best_score = s
            best_r = r

    # Fine search around best
    step = total_query / n_coarse
    fine_rots = np.linspace(best_r - 2 * step, best_r + 2 * step, 200)
    for r in fine_rots:
        s = monotonicity_score(int(r) % total_query)
        if s > best_score:
            best_score = s
            best_r = r

    rotate = round(best_r) % total_query
    return rotate, flip


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_dotplot(delta: DeltaFile, rotate: int = 0, flip: bool = False, ax=None):
    refs, queries, ref_offsets, query_offsets, total_ref, total_query = compute_layout(delta)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    else:
        fig = ax.figure

    rot_off = rotate % total_query if total_query > 0 else 0

    def to_display_y(cum_pos):
        """Map a cumulative query position to display position."""
        d = (cum_pos - rot_off) % total_query
        if flip:
            d = total_query - d
        return d

    # Draw alignments
    for sec in delta.sections:
        r_off = ref_offsets[sec.ref_id]
        q_off = query_offsets[sec.query_id]

        for aln in sec.alignments:
            is_forward = aln.query_start <= aln.query_end
            color = "#00BFFF" if is_forward else "#9933FF"

            x_start = r_off + aln.ref_start
            x_end = r_off + aln.ref_end
            cum_q_start = q_off + aln.query_start
            cum_q_end = q_off + aln.query_end

            if total_query == 0:
                continue

            disp_y_start = to_display_y(cum_q_start)
            disp_y_end = to_display_y(cum_q_end)

            # After rotation (pre-flip), check for wrapping
            rot_start = (cum_q_start - rot_off) % total_query
            rot_end = (cum_q_end - rot_off) % total_query
            wraps = (rot_start > rot_end) if is_forward else (rot_start < rot_end)

            if not wraps:
                ax.plot([x_start, x_end], [disp_y_start, disp_y_end],
                        color=color, linewidth=1, solid_capstyle="round")
            else:
                # Alignment wraps around — split into two segments
                if is_forward:
                    seg1 = float(total_query - rot_start)
                    seg2 = float(rot_end)
                else:
                    seg1 = float(rot_start)
                    seg2 = float(total_query - rot_end)
                fraction = seg1 / (seg1 + seg2) if (seg1 + seg2) > 0 else 0.5
                split_x = x_start + fraction * (x_end - x_start)

                # Forward alignments cross the total_query boundary in
                # rotated space; reverse alignments cross the 0 boundary.
                # Flip swaps which display edge each maps to.
                if is_forward != flip:
                    bnd1, bnd2 = total_query, 0
                else:
                    bnd1, bnd2 = 0, total_query

                ax.plot([x_start, split_x], [disp_y_start, bnd1],
                        color=color, linewidth=1, solid_capstyle="round")
                ax.plot([split_x, x_end], [bnd2, disp_y_end],
                        color=color, linewidth=1, solid_capstyle="round")

    # X-axis: auto-scaled numeric ticks
    def format_bp(val, _pos):
        return _format_bp(val)

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_bp))
    ax.tick_params(axis="x", labelsize=9)

    # Y-axis: place ticks at original (pre-rotation) coordinate values,
    # positioned at their display locations
    _nice_interval = _compute_nice_interval(total_query)
    orig_ticks = _generate_ticks(0, total_query, _nice_interval)
    # When rotating or flipping, drop tick at total_query to avoid overlap with 0
    if rot_off > 0 or flip:
        orig_ticks = [t for t in orig_ticks if round(t) != total_query]

    disp_positions = [to_display_y(round(t)) for t in orig_ticks]
    tick_labels = [_format_bp(t) for t in orig_ticks]

    ax.set_yticks(disp_positions)
    ax.set_yticklabels(tick_labels, fontsize=9)

    # Axis limits and labels
    ax.set_xlim(0, total_ref)
    ax.set_ylim(0, total_query)

    ref_label = delta.sections[0].ref_id if delta.sections else "Reference"
    ax.set_xlabel(ref_label, fontsize=11)

    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    return fig, ax


def _format_bp(val):
    if abs(val) >= 1_000_000:
        return f"{val / 1_000_000:.1f}M"
    elif abs(val) >= 1_000:
        return f"{val / 1_000:.1f}k"
    else:
        return f"{int(val)}"


def _compute_nice_interval(total):
    if total <= 0:
        return 1
    rough = total / 8
    magnitude = 10 ** int(np.floor(np.log10(rough)))
    normalized = rough / magnitude
    if normalized <= 1:
        return magnitude
    elif normalized <= 2:
        return 2 * magnitude
    elif normalized <= 5:
        return 5 * magnitude
    else:
        return 10 * magnitude


def _generate_ticks(start, end, interval):
    import math
    first = math.ceil(start / interval) * interval
    ticks = []
    t = first
    while t <= end + interval * 0.001:
        ticks.append(float(t))
        t += interval
    return ticks


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Dot plot visualisation for nucmer .delta files")
    parser.add_argument("delta", help="Path to .delta file")
    parser.add_argument("--rotate", type=int, default=0,
                        help="Circular rotation offset in bp for query (y-axis)")
    parser.add_argument("--flip", action="store_true",
                        help="Flip the query axis (reverse complement)")
    parser.add_argument("--auto", action="store_true",
                        help="Auto-detect optimal flip and rotation (overrides --flip and --rotate)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file (e.g. plot.png). If omitted, shows interactive window.")
    args = parser.parse_args()

    delta = parse_delta(args.delta)

    if args.auto:
        rotate, flip = auto_orient(delta)
        print(f"Auto-detected: rotate={rotate}, flip={flip}")
    else:
        rotate, flip = args.rotate, args.flip

    fig, ax = plot_dotplot(delta, rotate=rotate, flip=flip)
    plt.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
