#!/usr/bin/env python3
"""
got-dot-plot: An interactive DNA dot plot tool.

Compares two FASTA files using k-mer hashing and produces a self-contained
interactive HTML dot plot.

Usage:
    python dotplot.py x.fasta y.fasta
    python dotplot.py -k 15 -m 20 -o myplot.html x.fasta y.fasta
"""

import argparse
import logging
import sys
from pathlib import Path

import plotly.graph_objects as go
from Bio import SeqIO

logger = logging.getLogger(__name__)

# Byte-level complement table (upper-case only; N stays N)
_COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def parse_fasta(filepath):
    """Parse a FASTA file and return a list of (id, sequence_bytes) tuples."""
    path = Path(filepath)
    if not path.exists():
        logger.error("File not found: %s", filepath)
        sys.exit(1)

    records = [
        (r.id, bytes(str(r.seq).upper(), "ascii"))
        for r in SeqIO.parse(filepath, "fasta")
    ]
    if not records:
        logger.error("No sequences found in %s", filepath)
        sys.exit(1)
    return records


def concat_sequences(records, gap_size=100):
    """
    Concatenate sequences with N gaps so that no k-mer spans a boundary.

    Returns:
        combined   – bytes containing all sequences joined by gap_size N's
        boundaries – list of (name, start, end) where end is exclusive
    """
    gap = b"N" * gap_size
    buf = bytearray()
    boundaries = []

    for i, (name, seq) in enumerate(records):
        start = len(buf)
        buf.extend(seq)
        boundaries.append((name, start, len(buf)))
        if i < len(records) - 1:
            buf.extend(gap)

    return bytes(buf), boundaries


def reverse_complement(seq):
    """Return the reverse complement of a DNA bytes sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def build_kmer_index(seq, k):
    """
    Index every k-mer in seq.

    Returns a dict mapping kmer (bytes) -> list of start positions.
    K-mers containing 'N' are skipped.
    """
    idx = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if b"N" in kmer:
            continue
        if kmer in idx:
            idx[kmer].append(i)
        else:
            idx[kmer] = [i]
    return idx


def find_matches(x_seq, y_seq, k=10, min_length=10):
    """
    Find all exact-match diagonals between x_seq and y_seq.

    Algorithm:
      1. Build a k-mer index for both strands of x_seq.
      2. Walk every k-mer in y_seq and look it up in both indices.
      3. If the base before a hit also matches, this k-mer is already part
         of a longer match recorded at the previous position — skip it.
      4. Otherwise extend the match as far as bases agree and record it.

    Returns:
        forward_matches – list of (x_start, y_start, length)
        reverse_matches – list of (x_fwd_start, y_start, length)
            x_fwd_start is in forward-X coordinates; the match runs from
            (x_fwd_start, y_start) toward lower X and higher Y (anti-diagonal).
    """
    n_x = len(x_seq)
    n_y = len(y_seq)
    rc_x = reverse_complement(x_seq)
    rc_len = len(rc_x)

    logger.info("Indexing X forward strand...")
    fwd_idx = build_kmer_index(x_seq, k)
    logger.info("Indexing X reverse complement...")
    rev_idx = build_kmer_index(rc_x, k)

    forward_matches = []
    reverse_matches = []

    logger.info("Scanning Y sequence...")
    for i in range(n_y - k + 1):
        if i % 500_000 == 0 and i > 0:
            logger.debug("  %d%%  (%s / %s bp)", 100 * i // n_y, f"{i:,}", f"{n_y:,}")

        kmer = y_seq[i : i + k]
        if b"N" in kmer:
            continue

        # ---- Forward matches ----
        for x_pos in fwd_idx.get(kmer, ()):
            # If the preceding base also matched, this match was already
            # reported one position earlier — skip to avoid duplicates.
            if x_pos > 0 and i > 0 and x_seq[x_pos - 1] == y_seq[i - 1]:
                continue
            length = k
            while (
                x_pos + length < n_x
                and i + length < n_y
                and x_seq[x_pos + length] != ord("N")
                and y_seq[i + length] != ord("N")
                and x_seq[x_pos + length] == y_seq[i + length]
            ):
                length += 1
            if length >= min_length:
                forward_matches.append((x_pos, i, length))

        # ---- Reverse-complement matches ----
        for rc_pos in rev_idx.get(kmer, ()):
            if rc_pos > 0 and i > 0 and rc_x[rc_pos - 1] == y_seq[i - 1]:
                continue
            length = k
            while (
                rc_pos + length < rc_len
                and i + length < n_y
                and rc_x[rc_pos + length] != ord("N")
                and y_seq[i + length] != ord("N")
                and rc_x[rc_pos + length] == y_seq[i + length]
            ):
                length += 1
            if length >= min_length:
                # Map RC position → forward X coordinate.
                # The match runs anti-diagonally from (x_fwd, y_start)
                # toward lower X and higher Y.
                x_fwd = n_x - 1 - rc_pos
                reverse_matches.append((x_fwd, i, length))

    return forward_matches, reverse_matches


def _build_trace_coords(matches, forward=True):
    """
    Convert a match list to Plotly scatter coordinates.

    Segments are separated by None so they render as independent line
    segments within a single trace (far more efficient than one trace
    per match).

    Forward: (xStart, yStart) → (xStart+L, yStart+L)
    Reverse: (xFwd,   yStart) → (xFwd-L,   yStart+L)   [anti-diagonal]
    """
    xs, ys = [], []
    for x_start, y_start, length in matches:
        if forward:
            xs += [x_start, x_start + length, None]
            ys += [y_start, y_start + length, None]
        else:
            xs += [x_start, x_start - length, None]
            ys += [y_start, y_start + length, None]
    return xs, ys


def make_figure(x_boundaries, y_boundaries, fwd_matches, rev_matches, x_label, y_label):
    """Build and return the Plotly Figure."""
    x_total = x_boundaries[-1][2]
    y_total = y_boundaries[-1][2]

    traces = []

    if fwd_matches:
        xs, ys = _build_trace_coords(fwd_matches, forward=True)
        traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color="steelblue", width=1),
                name="Forward",
            )
        )

    if rev_matches:
        xs, ys = _build_trace_coords(rev_matches, forward=False)
        traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color="crimson", width=1),
                name="Reverse complement",
            )
        )

    shapes = []
    annotations = []

    # Vertical dashed lines between X sequences
    for _name, _start, end in x_boundaries[:-1]:
        shapes.append(
            dict(
                type="line",
                x0=end,
                x1=end,
                y0=0,
                y1=y_total,
                line=dict(color="#cccccc", width=1, dash="dash"),
                layer="below",
            )
        )

    # Horizontal dashed lines between Y sequences
    for _name, _start, end in y_boundaries[:-1]:
        shapes.append(
            dict(
                type="line",
                x0=0,
                x1=x_total,
                y0=end,
                y1=end,
                line=dict(color="#cccccc", width=1, dash="dash"),
                layer="below",
            )
        )

    # Sequence name labels along X axis (above plot area)
    for name, start, end in x_boundaries:
        annotations.append(
            dict(
                x=(start + end) / 2,
                y=y_total,
                xref="x",
                yref="y",
                text=f"<b>{name}</b>",
                showarrow=False,
                xanchor="center",
                yanchor="bottom",
                font=dict(size=9),
            )
        )

    # Sequence name labels along Y axis (left of plot area)
    for name, start, end in y_boundaries:
        annotations.append(
            dict(
                x=0,
                y=(start + end) / 2,
                xref="x",
                yref="y",
                text=f"<b>{name}</b>",
                showarrow=False,
                xanchor="right",
                yanchor="middle",
                font=dict(size=9),
            )
        )

    fig = go.Figure(data=traces)
    fig.update_layout(
        xaxis=dict(
            title=x_label,
            range=[0, x_total],
            minallowed=0,
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            title=y_label,
            range=[0, y_total],
            minallowed=0,
            showgrid=False,
            zeroline=False,
            scaleanchor="x",  # equal aspect ratio: diagonals appear at 45°
            scaleratio=1,
        ),
        shapes=shapes,
        annotations=annotations,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        margin=dict(l=120, r=20, t=60, b=60),
    )

    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive DNA dot plot from two FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("x_fasta", help="FASTA file for the X axis")
    parser.add_argument("y_fasta", help="FASTA file for the Y axis")
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=10, help="K-mer size for seeding matches"
    )
    parser.add_argument(
        "-m",
        "--min-match",
        type=int,
        default=100,
        help="Minimum match length to report",
    )
    parser.add_argument(
        "-o", "--output", default="dotplot.html", help="Output HTML file"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show debug-level progress messages",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    logger.info("Parsing %s...", args.x_fasta)
    x_records = parse_fasta(args.x_fasta)
    logger.info(
        "  %d sequence(s), %s bp total",
        len(x_records),
        f"{sum(len(s) for _, s in x_records):,}",
    )

    logger.info("Parsing %s...", args.y_fasta)
    y_records = parse_fasta(args.y_fasta)
    logger.info(
        "  %d sequence(s), %s bp total",
        len(y_records),
        f"{sum(len(s) for _, s in y_records):,}",
    )

    x_seq, x_bounds = concat_sequences(x_records)
    y_seq, y_bounds = concat_sequences(y_records)

    logger.info(
        "Finding matches (k=%d, min_length=%d)...", args.kmer_size, args.min_match
    )
    fwd, rev = find_matches(x_seq, y_seq, k=args.kmer_size, min_length=args.min_match)
    logger.info(
        "  %s forward match(es), %s reverse match(es)", f"{len(fwd):,}", f"{len(rev):,}"
    )

    logger.info("Building plot...")
    fig = make_figure(
        x_bounds, y_bounds, fwd, rev, Path(args.x_fasta).stem, Path(args.y_fasta).stem
    )

    fig.write_html(args.output)
    logger.info("Done!  Open %s in your browser.", args.output)


if __name__ == "__main__":
    main()
