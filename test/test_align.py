# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    nw.align(seq1, seq2)

    assert nw._align_matrix is not None
    assert nw._gapA_matrix is not None
    assert nw._gapB_matrix is not None

    expected_align_matrix = np.array([[0, -10, -11, -12],[-10, 5, -5, -6],[-11, -5, -1, -2],[-12, -6, 10, 0],[-13, -7, 0, 15]])
    expected_gapA_matrix = np.array([[0, -10, -20, -30],[-10, -15, -25, -35],[-20, -25, -30, -40],[-30, -35, -45, -50],[-40, -45, -55, -60]])
    expected_gapB_matrix = np.array([[0, -10, -11, -12],[-10, -15, -16, -17],[-11, -16, -17, -18],[-12, -17, -18, -19],[-13, -18, -19, -20]])

    assert (nw._align_matrix == expected_align_matrix).all(), "Alignment matrix does not match expected values."
    assert (nw._gapA_matrix == expected_gapA_matrix).all(), "GapA matrix does not match expected values."
    assert (nw._gapB_matrix == expected_gapB_matrix).all(), "GapB matrix does not match expected values."

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, aligned_seq3, aligned_seq4 = nw.align(seq3, seq4)

    # Expected results based on given example
    expected_score = 17
    expected_aligned_seq3 = "MAVHQLIRRP"
    expected_aligned_seq4 = "M---QLIRHP"

    # Assertions to validate alignment
    assert score == expected_score, f"Expected alignment score {expected_score}, but got {score}."
    assert aligned_seq3 == expected_aligned_seq3, f"Expected aligned sequence {expected_aligned_seq3}, but got {aligned_seq3}."
    assert aligned_seq4 == expected_aligned_seq4, f"Expected aligned sequence {expected_aligned_seq4}, but got {aligned_seq4}."
