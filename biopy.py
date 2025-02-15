from Bio.Align import substitution_matrices
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# load sub matrix BLOSUM62
matrix = substitution_matrices.load("BLOSUM62")
# sequence to be aligned
seq1 = "MAVHQLIRRP"
seq2 = "MQLIRHP"

# gap penalty
gap_open = -10
gap_extend = -1


# global alignment
alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

# alignment result
for alignment in alignments:
    print(format_alignment(*alignment))
