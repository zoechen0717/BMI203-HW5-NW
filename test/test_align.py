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

    # Assert the initialized matrices are correct
    assert nw._align_matrix is not None
    assert nw._gapA_matrix is not None
    assert nw._gapB_matrix is not None

    # Expected align matrix generated from test.py
    expected_align_matrix = np.array([
        [[  0. -10. -11. -12.],
        [-10.   5.  -5.  -6.],
        [-11.  -5.   4.  -6.],
        [-12.  -6.   0.   5.],
        [-13.  -7.  -5.   5.]]
    expected_gapA_matrix = np.array([
        [[-inf -inf -inf -inf],
        [-10. -20. -21. -22.],
        [-11.  -5. -15. -16.],
        [-12.  -6.  -6. -16.],
        [-13.  -7.  -7.  -5.]]

    expected_gapB_matrix = np.array([
        Gap B Matrix:
        [[-inf -10. -11. -12.],
        [-inf -20.  -5.  -6.],
        [-inf -21. -15.  -6.],
        [-inf -22. -16. -10.],
        [-inf -23. -17. -15.]]
    # Assert the output matrices
    assert (nw._align_matrix == expected_align_matrix).all(), "Alignment matrix does not match."
    assert (nw._gapA_matrix == expected_gapA_matrix).all(), "GapA matrix does not match."
    assert (nw._gapB_matrix == expected_gapB_matrix).all(), "GapB matrix does not match."

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
    # Result from biopython pairwise2.align.globalds equal 18 not 17
    # Check with the biopy.py
    expected_score = 18
    # Expected socre and backtracked seq are generated using biopython pairwise2.align.globalds
    expected_aligned_seq3 = "MAVHQLIRRP"
    expected_aligned_seq4 = "M---QLIRHP"

    # Assertions to validate alignment
    assert score == expected_score, f"Incorrect alignment score."
    assert aligned_seq3 == expected_aligned_seq3, f"Incorrect aligned sequence."
    assert aligned_seq4 == expected_aligned_seq4, f"Incorrect aligned sequence."
