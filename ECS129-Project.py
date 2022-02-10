print("Hello")
input = "TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT"

"""
Reverse Compliment DNA
Created by Madison Hypes
def reverse_compliment(seq):
"""


def reverse_compliment(seq):

    return seq


"""
This is some extra stuff from the tutorial. I can recode it.
inputfile = "DNA_sequence_original.txt"
f = open(inputfile, "r")
seq = f.read()

seq = seq.replace("\n", "")
seq = seq.replace("\r", "")
"""

"""
Created by arorakashish0911 on geekforgeeks.org tutorials
https://www.geeksforgeeks.org/dna-protein-python-3/
def translate(seq)
"""


def translate(seq):
    protein = ""
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        protein += table[codon]
    return protein


"""
Created by arorakashish0911 on geekforgeeks.org tutorials
https://www.geeksforgeeks.org/dna-protein-python-3/
def read_seq(inputfile):
"""


def read_seq(inputfile):
    with open(inputfile, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


"""
Extra tutorial stuff
prt = read_seq("amino_acid_sequence_original.txt")
dna = read_seq("DNA_sequence_original.txt")


p = translate(dna[20:935])
p == prt
"""
