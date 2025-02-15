# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # Prepare dic for looping
    species_sequences = {
        gg_header: gg_seq,
        mm_header: mm_seq,
        br_header: br_seq,
        tt_header: tt_seq}
    # Initialize Needleman-Wunsch with BLOSUM62 matrix and gap penalties
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # Looping perform alignments with human sequences and store scores
    scores = {}
    for species, seq in species_sequences.items():
        score, _, _ = nw.align(hs_seq, seq)
        scores[species] = score

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # Sort species by similarity (highest score first)
    sorted_species = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    # Print sorted species and scores
    for species, score in sorted_species:
        print(f"{species}: Alignment score = {score:.2f}")


if __name__ == "__main__":
    main()
