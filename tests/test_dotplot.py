"""Tests for dotplot.py"""

from pathlib import Path

import plotly.graph_objects as go
import pytest

from dotplot import (
    build_kmer_index,
    concat_sequences,
    find_matches,
    make_figure,
    parse_fasta,
    reverse_complement,
)

# ---------------------------------------------------------------------------
# Helper: a simple RC implementation independent of the code under test,
# used only to set up test fixtures.
# ---------------------------------------------------------------------------


def _rc(seq: bytes) -> bytes:
    comp = bytes.maketrans(b"ACGTN", b"TGCAN")
    return seq.translate(comp)[::-1]


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------


class TestReverseComplement:
    def test_single_bases(self):
        assert reverse_complement(b"A") == b"T"
        assert reverse_complement(b"T") == b"A"
        assert reverse_complement(b"C") == b"G"
        assert reverse_complement(b"G") == b"C"

    def test_n_preserved(self):
        assert reverse_complement(b"N") == b"N"

    def test_known_sequences(self):
        assert reverse_complement(b"ACGT") == b"ACGT"  # palindrome
        assert reverse_complement(b"AAAA") == b"TTTT"
        assert reverse_complement(b"AACCGG") == b"CCGGTT"

    def test_involution(self):
        """Applying RC twice returns the original sequence."""
        seq = b"ACGATCGATCACGATCGATC"
        assert reverse_complement(reverse_complement(seq)) == seq


# ---------------------------------------------------------------------------
# concat_sequences
# ---------------------------------------------------------------------------


class TestConcatSequences:
    def test_single_sequence_no_gap(self):
        records = [("seq1", b"ACGT")]
        combined, boundaries = concat_sequences(records, gap_size=10)
        assert combined == b"ACGT"
        assert boundaries == [("seq1", 0, 4)]

    def test_two_sequences_with_gap(self):
        records = [("a", b"AAAA"), ("b", b"CCCC")]
        combined, boundaries = concat_sequences(records, gap_size=5)
        assert combined == b"AAAA" + b"N" * 5 + b"CCCC"
        assert boundaries == [("a", 0, 4), ("b", 9, 13)]

    def test_three_sequences(self):
        records = [("a", b"AA"), ("b", b"CC"), ("c", b"GG")]
        combined, boundaries = concat_sequences(records, gap_size=2)
        assert combined == b"AANNCCNNGG"
        assert boundaries == [("a", 0, 2), ("b", 4, 6), ("c", 8, 10)]

    def test_boundary_end_is_exclusive(self):
        records = [("s", b"ACGT")]
        _, boundaries = concat_sequences(records)
        _name, start, end = boundaries[0]
        assert end - start == 4


# ---------------------------------------------------------------------------
# build_kmer_index
# ---------------------------------------------------------------------------


class TestBuildKmerIndex:
    def test_basic_indexing(self):
        idx = build_kmer_index(b"ACGTACGT", k=4)
        assert b"ACGT" in idx
        assert idx[b"ACGT"] == [0, 4]

    def test_n_containing_kmers_excluded(self):
        # ACNGT has 3-mers: ACN, CNG, NGT — all contain N
        idx = build_kmer_index(b"ACNGT", k=3)
        assert len(idx) == 0

    def test_n_excluded_only_at_boundary(self):
        # ACNACGT has 3-mers: ACN(skip), CNA(skip), NAC(skip), ACG, CGT
        idx = build_kmer_index(b"ACNACGT", k=3)
        assert b"ACG" in idx
        assert b"CGT" in idx
        assert idx[b"ACG"] == [3]
        assert idx[b"CGT"] == [4]

    def test_unique_kmer_single_position(self):
        idx = build_kmer_index(b"ACGATCG", k=4)
        assert idx[b"ACGA"] == [0]

    def test_duplicate_kmer_multiple_positions(self):
        # 9-char sequence so ACGT appears at 0 and 4, CGTA at 1 and 5
        idx = build_kmer_index(b"ACGTACGTA", k=4)
        assert idx[b"ACGT"] == [0, 4]
        assert idx[b"CGTA"] == [1, 5]


# ---------------------------------------------------------------------------
# find_matches
# ---------------------------------------------------------------------------


class TestFindMatches:
    # 20 bp motifs used across several tests.
    # Verified non-palindromic (RC ≠ self) and non-repetitive (no repeated 10-mers).
    FWD_MOTIF = b"AACGTTAGCGATCATGCCAT"
    RC_SRC = b"GCTATCAGTTACGGATCTAG"

    def test_forward_match_found(self):
        """A shared motif produces exactly one forward match at the expected coords."""
        x = b"GGGGGGGGGG" + self.FWD_MOTIF + b"CCCCCCCCCC"
        y = b"TTTTTTTTTT" + self.FWD_MOTIF + b"AAAAAAAAAA"
        fwd, rev = find_matches(x, y, k=10, min_length=10)
        assert (10, 10, 20) in fwd

    def test_no_matches_for_different_sequences(self):
        x = b"A" * 30
        y = b"C" * 30
        fwd, rev = find_matches(x, y, k=10, min_length=10)
        assert fwd == []
        assert rev == []

    def test_reverse_complement_match_found(self):
        """RC of a motif in X appearing in Y produces a reverse match."""
        rc_src = _rc(self.RC_SRC)
        x = b"GGGGGGGGGG" + self.RC_SRC + b"CCCCCCCCCC"  # 40 bp
        y = b"TTTTTTTTTT" + rc_src + b"AAAAAAAAAA"  # 40 bp
        fwd, rev = find_matches(x, y, k=10, min_length=10)
        # rc_pos = 10, i = 10, length = 20
        # x_fwd = len(x) - 1 - rc_pos = 40 - 1 - 10 = 29
        assert (29, 10, 20) in rev
        assert fwd == []

    def test_min_length_filters_short_matches(self):
        """Matches shorter than min_length are not reported."""
        motif = b"AACGTTAGCGATCAT"  # 15 bp
        x = b"GGGGGGGGGG" + motif + b"CCCCCCCCCC"
        y = b"TTTTTTTTTT" + motif + b"AAAAAAAAAA"
        fwd_included, _ = find_matches(x, y, k=10, min_length=10)
        fwd_excluded, _ = find_matches(x, y, k=10, min_length=20)
        assert any(xs == 10 for xs, _, _ in fwd_included)
        assert not any(xs == 10 for xs, _, _ in fwd_excluded)

    def test_no_duplicate_matches(self):
        """Each distinct match region is reported exactly once."""
        x = b"GGGGGGGGGG" + self.FWD_MOTIF + b"CCCCCCCCCC"
        y = b"TTTTTTTTTT" + self.FWD_MOTIF + b"AAAAAAAAAA"
        fwd, _ = find_matches(x, y, k=10, min_length=10)
        positions = [(xs, ys) for xs, ys, _ in fwd]
        assert len(positions) == len(set(positions))

    def test_identical_sequences(self):
        """Identical sequences produce one forward match spanning the full length."""
        seq = b"ACGTACGTTGCAAGCTATGC"  # 20 bp, no repeated 10-mers
        fwd, _rev = find_matches(seq, seq, k=10, min_length=10)
        assert fwd == [(0, 0, 20)]


# ---------------------------------------------------------------------------
# make_figure
# ---------------------------------------------------------------------------


class TestMakeFigure:
    X_BOUNDS = [("chrX", 0, 1000)]
    Y_BOUNDS = [("chrY", 0, 1000)]

    def test_returns_plotly_figure(self):
        fig = make_figure(self.X_BOUNDS, self.Y_BOUNDS, [], [], "x", "y")
        assert isinstance(fig, go.Figure)

    def test_no_traces_for_empty_matches(self):
        fig = make_figure(self.X_BOUNDS, self.Y_BOUNDS, [], [], "x", "y")
        assert len(fig.data) == 0

    def test_forward_trace_present(self):
        fig = make_figure(self.X_BOUNDS, self.Y_BOUNDS, [(0, 0, 100)], [], "x", "y")
        assert len(fig.data) == 1
        assert fig.data[0].name == "Forward"

    def test_reverse_trace_present(self):
        fig = make_figure(self.X_BOUNDS, self.Y_BOUNDS, [], [(500, 0, 100)], "x", "y")
        assert len(fig.data) == 1
        assert fig.data[0].name == "Reverse complement"

    def test_both_traces_present(self):
        fig = make_figure(
            self.X_BOUNDS, self.Y_BOUNDS, [(0, 0, 100)], [(500, 0, 100)], "x", "y"
        )
        assert len(fig.data) == 2

    def test_boundary_shapes_for_multi_sequence(self):
        """One vertical and one horizontal dashed line for two-sequence inputs."""
        x_bounds = [("a", 0, 500), ("b", 600, 1000)]
        y_bounds = [("c", 0, 400), ("d", 500, 1000)]
        fig = make_figure(x_bounds, y_bounds, [], [], "x", "y")
        assert len(fig.layout.shapes) == 2


# ---------------------------------------------------------------------------
# parse_fasta
# ---------------------------------------------------------------------------


class TestParseFasta:
    def test_single_sequence(self, tmp_path):
        f = tmp_path / "test.fasta"
        f.write_text(">seq1\nACGTACGT\n")
        records = parse_fasta(f)
        assert records == [("seq1", b"ACGTACGT")]

    def test_multiple_sequences(self, tmp_path):
        f = tmp_path / "test.fasta"
        f.write_text(">a\nAAAA\n>b\nCCCC\n")
        records = parse_fasta(f)
        assert len(records) == 2
        assert records[0] == ("a", b"AAAA")
        assert records[1] == ("b", b"CCCC")

    def test_lowercase_is_uppercased(self, tmp_path):
        f = tmp_path / "test.fasta"
        f.write_text(">seq1\nacgtacgt\n")
        records = parse_fasta(f)
        assert records[0][1] == b"ACGTACGT"

    def test_missing_file_exits(self):
        with pytest.raises(SystemExit):
            parse_fasta("/nonexistent/path/file.fasta")

    def test_empty_file_exits(self, tmp_path):
        f = tmp_path / "empty.fasta"
        f.write_text("")
        with pytest.raises(SystemExit):
            parse_fasta(f)


# ---------------------------------------------------------------------------
# Integration: toy 1000 bp genomes
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def toy_matches():
    x_rec = parse_fasta(DATA_DIR / "seq_x.fasta")
    y_rec = parse_fasta(DATA_DIR / "seq_y.fasta")
    x_seq, _ = concat_sequences(x_rec)
    y_seq, _ = concat_sequences(y_rec)
    fwd, rev = find_matches(x_seq, y_seq, k=10, min_length=50)
    return fwd, rev


class TestToyGenomes:
    """Integration tests against the 1000 bp toy FASTA files in tests/data/.

    The sequences were generated with known structure:
      - Forward match:  chrX[100:300] == chrY[100:300]  → (100, 100, 200)
      - Reverse match:  RC(chrX[300:500]) == chrY[400:600] → (499, 400, 200)

    min_length=50 is used to rule out any accidental short matches from the
    random flanking regions.
    """

    def test_forward_match_found(self, toy_matches):
        fwd, _rev = toy_matches
        assert (100, 100, 200) in fwd

    def test_reverse_match_found(self, toy_matches):
        _fwd, rev = toy_matches
        assert (499, 400, 200) in rev

    def test_html_output_written(self, toy_matches, tmp_path):
        x_rec = parse_fasta(DATA_DIR / "seq_x.fasta")
        y_rec = parse_fasta(DATA_DIR / "seq_y.fasta")
        _, x_bounds = concat_sequences(x_rec)
        _, y_bounds = concat_sequences(y_rec)
        fwd, rev = toy_matches
        fig = make_figure(x_bounds, y_bounds, fwd, rev, "seq_x", "seq_y")
        out = tmp_path / "dotplot.html"
        fig.write_html(str(out))
        assert out.exists()
        assert out.stat().st_size > 0
