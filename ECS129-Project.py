import os

os.chdir(r'C:\Users\Madison Hypes\Desktop\FOLDER\Code\ECS129Project')

# Created by arorakashish0911 on geekforgeeks.org tutorials
# https: // www.geeksforgeeks.org/dna-protein-python-3/

# Python program to convert
# altered DNA to protein

# inputfile = open("DNAfastaExample.txt", "r")
# f = inputfile
# seq = f.read()
# seq = seq.replace("\n", "")
# seq = seq.replace("\r", "")
# f.close()
# print(seq)

# Function created by Madison Hypes


def readFasta():
    with open("DNAfastaExample.txt", "r") as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                pass
            else:
                line = line.replace("\n", "")
                line = line.replace("\r", "")
                line = line.replace(" ", "")
                seq += line.strip()
        f.close()
    return seq


seq = readFasta()
print(seq)


def translate(seq):

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
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


def read_seq(inputfile):
    with open(inputfile, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


# prt = read_seq("amino_acid_sequence_original.txt")
# dna = read_seq(inputfile)
# print(dna)
# p = translate(dna[20:935])
# p == prt
