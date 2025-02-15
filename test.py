import numpy as np
from align import NeedlemanWunsch, read_fasta

def test_local():
    """
    Runs a local test for the Needleman-Wunsch alignment.
    """
    # Load test sequences
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    # Initialize Needleman-Wunsch with BLOSUM62 and gap penalties
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # Perform alignment
    score, aligned_seq1, aligned_seq2 = nw.align(seq1, seq2)

    # Print results
    print("\n--- Needleman-Wunsch Local Test ---")
    print(f"Input Sequence 1: {seq1}")
    print(f"Input Sequence 2: {seq2}")
    print("\nAlignment Matrices:")
    print("\nAlignment Score Matrix:")
    print(nw._align_matrix)
    print("\nGap A Matrix:")
    print(nw._gapA_matrix)
    print("\nGap B Matrix:")
    print(nw._gapB_matrix)
    print("\nBacktracking Matrix:")
    print(nw._back)

    print("\nFinal Alignment:")
    print(f"Score: {score}")
    print(f"Aligned Sequence 1: {aligned_seq1}")
    print(f"Aligned Sequence 2: {aligned_seq2}")

# Run the test
if __name__ == "__main__":
    test_local()
