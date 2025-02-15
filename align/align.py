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
        TODO

        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm

        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA

        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # The align function and backtrace functions are inspired by https://github.com/ahishsujay/Sequence_Alignment/blob/master/nwAlign.py
        # ChatGPT was used for all the debugging steps ---> all the print functions

        # Reset alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Reset alignment score in case method is called more than once
        self.alignment_score = 0

        self._seqA = seqA
        self._seqB = seqB
        m, n = len(seqA), len(seqB)

        # Initialize alignment and gap matrices
        self._align_matrix = np.zeros((m + 1, n + 1))
        self._gapA_matrix = np.full((m + 1, n + 1), float('-inf'))
        self._gapB_matrix = np.full((m + 1, n + 1), float('-inf'))
        #self._gapA_matrix = np.zeros((m + 1, n + 1))
        #self._gapB_matrix = np.zeros((m + 1, n + 1))
        self._back = np.zeros((m + 1, n + 1), dtype=int)

        #print("\n*** INITIALIZING MATRICES ***\n")

        # Initialize gap penalties in the first row and column
        for i in range(1, m + 1):
            self._align_matrix[i][0] = self.gap_open + (i - 1) * self.gap_extend
            self._back[i][0] = 1  # Indicates a gap in seqB
            #print(f"Init align_matrix[{i}][0]: {self._align_matrix[i][0]} (Gap in seqB)")

        for j in range(1, n + 1):
            self._align_matrix[0][j] = self.gap_open + (j - 1) * self.gap_extend
            self._back[0][j] = 2  # Indicates a gap in seqA
            #print(f"Init align_matrix[0][{j}]: {self._align_matrix[0][j]} (Gap in seqA)")

        # Initialize first column (seqB gaps)
        for i in range(1, m + 1):
            self._gapA_matrix[i][0] = self.gap_open + (i - 1) * self.gap_extend

        # Initialize first row (seqA gaps)
        for j in range(1, n + 1):
            self._gapB_matrix[0][j] = self.gap_open + (j - 1) * self.gap_extend

        self._align_matrix[0][0] = 0  # Ensure top-left cell is 0
        #print("\nAlignment Matrix Initialized:\n", self._align_matrix)

        # Fill alignment and backtracking matrices
        #print("\n*** FILLING MATRICES ***\n")

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Match/mismatch case
                match_score = self.sub_dict.get((seqA[i-1], seqB[j-1]), -1)
                match = self._align_matrix[i-1][j-1] + match_score
                # Handling gaps properly
                gapA = max(self._align_matrix[i-1][j] + self.gap_open, self._gapA_matrix[i-1][j] + self.gap_extend)
                gapB = max(self._align_matrix[i][j-1] + self.gap_open, self._gapB_matrix[i][j-1] + self.gap_extend)
                # Fill the alignment matrix with the best score
                self._align_matrix[i][j] = max(match, gapA, gapB)
                # Correctly update gap matrices
                self._gapA_matrix[i][j] = gapA
                self._gapB_matrix[i][j] = gapB
                # print score calculations for debugging
                #print(f"Cell ({i}, {j}): Match={match}, GapA={gapA}, GapB={gapB} -> Max={self._align_matrix[i][j]}")

                # Backtracking logic
                if self._align_matrix[i][j] == match:
                    self._back[i][j] = 0  # Match/mismatch
                    #print(f"Backtracking: ({i}, {j}) → Match")
                elif self._align_matrix[i][j] == gapA:
                    self._back[i][j] = 1  # Gap in seqB
                    #print(f"Backtracking: ({i}, {j}) → Gap in seqB")
                else:
                    self._back[i][j] = 2  # Gap in seqA
                    #print(f"Backtracking: ({i}, {j}) → Gap in seqA")

        #print("\nFinal Alignment Matrix:\n", self._align_matrix)
        #print("\nBacktracking Matrix:\n", self._back)

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO

        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.

        Parameters:
        	None

        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # The align function and backtrace functions are inspired by https://github.com/ahishsujay/Sequence_Alignment/blob/master/nwAlign.py
        # initializing matrix need for backtracing
        m, n = len(self._seqA), len(self._seqB)
        i, j = m, n
        self.alignment_score = self._align_matrix[i][j]
        seqA_align, seqB_align = [], []

        while i > 0 or j > 0:
            # trace back for match
            if i > 0 and j > 0 and self._back[i][j] == 0:
                seqA_align.append(self._seqA[i-1])
                seqB_align.append(self._seqB[j-1])
                i -= 1
                j -= 1
            # trace back for gap in B
            elif i > 0 and self._back[i][j] == 1:
                seqA_align.append(self._seqA[i-1])
                seqB_align.append('-')
                i -= 1
            # trace back for gap in A
            else:
                seqA_align.append('-')
                seqB_align.append(self._seqB[j-1])
                j -= 1
        self.seqA_align = ''.join(reversed(seqA_align))
        self.seqB_align = ''.join(reversed(seqB_align))

        return self.alignment_score, self.seqA_align, self.seqB_align


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
