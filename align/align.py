# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/fimme of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm with affine gap penalties

        Parameters:
            seqA: str
                the first string to be aligned
            seqB: str
                the second string to be aligned with seqA

        Returns:
            (alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
                the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""
        self.alignment_score = 0
        self._seqA = seqA
        self._seqB = seqB

        lenA, lenB = len(seqA), len(seqB)

        # Initialize score matrices
        self._align_matrix = np.full((lenA+1, lenB+1), -np.inf)
        self._gapA_matrix = np.full((lenA+1, lenB+1), -np.inf)  # Gaps in A
        self._gapB_matrix = np.full((lenA+1, lenB+1), -np.inf)  # Gaps in B

        # Initialize traceback matrices
        self._back = np.zeros((lenA+1, lenB+1), dtype=object)
        self._back_A = np.zeros((lenA+1, lenB+1), dtype=object)
        self._back_B = np.zeros((lenA+1, lenB+1), dtype=object)

        # Base cases
        self._align_matrix[0][0] = 0

        # Initialize gapA matrix first row (gaps in A)
        for j in range(1, lenB+1):
            self._gapA_matrix[0][j] = self.gap_open + (j-1)*self.gap_extend
            self._back_A[0][j] = ('gapA', 0, j-1)

        # Initialize gapB matrix first column (gaps in B)
        for i in range(1, lenA+1):
            self._gapB_matrix[i][0] = self.gap_open + (i-1)*self.gap_extend
            self._back_B[i][0] = ('gapB', i-1, 0)

        # Fill matrices
        for i in range(1, lenA+1):
            for j in range(1, lenB+1):
                # Update match/mismatch matrix
                score = self.sub_dict.get((seqA[i-1], seqB[j-1]), 0)

                align_prev = self._align_matrix[i-1][j-1] + score
                gapA_prev = self._gapA_matrix[i-1][j-1] + score
                gapB_prev = self._gapB_matrix[i-1][j-1] + score
                max_val = max(align_prev, gapA_prev, gapB_prev)
                self._align_matrix[i][j] = max_val

                # Record traceback
                if max_val == align_prev:
                    self._back[i][j] = ('align', i-1, j-1)
                elif max_val == gapA_prev:
                    self._back[i][j] = ('gapA', i-1, j-1)
                else:
                    self._back[i][j] = ('gapB', i-1, j-1)

                # Update gapA matrix (gaps in A)
                extend = self._gapA_matrix[i][j-1] + self.gap_extend
                new_gap = self._align_matrix[i][j-1] + self.gap_open
                max_gapA = max(extend, new_gap)
                self._gapA_matrix[i][j] = max_gapA

                if max_gapA == extend:
                    self._back_A[i][j] = ('gapA', i, j-1)
                else:
                    self._back_A[i][j] = ('align', i, j-1)

                # Update gapB matrix (gaps in B)
                extend = self._gapB_matrix[i-1][j] + self.gap_extend
                new_gap = self._align_matrix[i-1][j] + self.gap_open
                max_gapB = max(extend, new_gap)
                self._gapB_matrix[i][j] = max_gapB

                if max_gapB == extend:
                    self._back_B[i][j] = ('gapB', i-1, j)
                else:
                    self._back_B[i][j] = ('align', i-1, j)

        # Determine final score
        self.alignment_score = max(self._align_matrix[lenA][lenB],
                                  self._gapA_matrix[lenA][lenB],
                                  self._gapB_matrix[lenA][lenB])

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        Trace back through the matrices to construct the alignment

        Returns:
            (alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
        """
        i, j = len(self._seqA), len(self._seqB)
        alignA, alignB = [], []

        # Determine starting matrix
        current = None
        max_score = max(self._align_matrix[i][j],
                        self._gapA_matrix[i][j],
                        self._gapB_matrix[i][j])

        if max_score == self._align_matrix[i][j]:
            current = ('align', i, j)
        elif max_score == self._gapA_matrix[i][j]:
            current = ('gapA', i, j)
        else:
            current = ('gapB', i, j)

        while i > 0 or j > 0:
            matrix_type, i, j = current

            if matrix_type == 'align':
                alignA.append(self._seqA[i-1])
                alignB.append(self._seqB[j-1])
                current = self._back[i][j]
                i -= 1
                j -= 1
            elif matrix_type == 'gapA':
                alignA.append('-')
                alignB.append(self._seqB[j-1])
                current = self._back_A[i][j]
                j -= 1
            elif matrix_type == 'gapB':
                alignA.append(self._seqA[i-1])
                alignB.append('-')
                current = self._back_B[i][j]
                i -= 1

        # Reverse alignments
        self.seqA_align = ''.join(reversed(alignA))
        self.seqB_align = ''.join(reversed(alignB))

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
