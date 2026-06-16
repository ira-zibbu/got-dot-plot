# got-dot-plot

An interactive DNA dot plot tool. Compares two FASTA files and produces a
self-contained HTML file you can open in any browser — zoom, pan, and export
without any server or extra software.

Built as a modern replacement for [Re-dot-able](https://github.com/s-andrews/redotable),
without the Java dependency headaches.

## What is a dot plot?

A dot plot is a way to visually compare two DNA sequences. Each axis represents
one sequence, and a diagonal line means that region is identical (or
reverse-complementary) between the two. You can immediately see:

- **Conserved regions** — long diagonal lines
- **Inversions** — anti-diagonal lines (reverse-complement matches)
- **Rearrangements** — diagonal lines that are offset or out of order
- **Duplications** — multiple parallel diagonals

## Installation

Requires Python ≥ 3.14. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv tool install --editable .
```

This makes `got-dot-plot` available as a command anywhere on your system.

## Usage

```bash
got-dot-plot reference.fasta assembly.fasta
```

This writes `dotplot.html` to the current directory. Open it in your browser.

```bash
# Custom output file
got-dot-plot reference.fasta assembly.fasta -o myplot.html

# Increase minimum match length for cleaner plots on large genomes
got-dot-plot reference.fasta assembly.fasta -m 500

# Show verbose progress
got-dot-plot reference.fasta assembly.fasta -v
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-k`, `--kmer-size` | `10` | K-mer size for seeding matches. Larger values are faster and more specific but miss divergent regions. |
| `-m`, `--min-match` | `100` | Minimum match length to display. Increase this to reduce HTML file size and improve browser performance. |
| `-o`, `--output` | `dotplot.html` | Output HTML file. |
| `-v`, `--verbose` | off | Show per-500 kbp progress messages during scanning. |

## How it works

Matching uses k-mer hashing — the same approach as
[MUMmer](https://mummer4.github.io/) dot plots:

1. Every k-mer in the X sequence (both forward and reverse-complement strands)
   is indexed into a hash table.
2. Each k-mer in the Y sequence is looked up in the index.
3. On a hit, the match is extended in both directions until a mismatch or N is
   reached.
4. Matches shorter than `--min-match` are discarded.

**Forward matches** (same orientation) appear as blue diagonal lines.
**Reverse-complement matches** (inversions) appear as red anti-diagonal lines.

Both FASTA files can contain multiple sequences. They are concatenated with
N gaps for display, with sequence names labelled on the axes and dashed lines
marking boundaries.

## Tuning for performance

| Genome size | Recommended `-m` |
|-------------|-----------------|
| < 100 kbp | 50 |
| 100 kbp – 10 Mbp | 100–500 |
| > 10 Mbp | 500+ |

Increasing `-k` also speeds up the scan by reducing spurious k-mer hits, at
the cost of missing matches that contain mismatches in the seed region.

## Development

```bash
# Install dependencies including dev tools
uv sync --group dev

# Run tests
uv run pytest
```
