import os
from tracemalloc import start

os.chdir(r'C:\Users\Madison Hypes\Desktop\FOLDER\Code\ECS129Project')


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


def dnaPol(seq):
    table = {"A": "T", "T": "A", "G": "C", "C": "G"}
    reverseComp = ""
    for i in range(0, len(seq), 1):
        bp = seq[i]
        reverseComp += table[bp]
    reverseComp = ''.join(reversed(reverseComp))
    return reverseComp


def rnaPol():
    frame1 = ""
    frame2 = ""
    frame3 = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon == "ATG":
            start = i
            for i in range(start, len(seq), 3):
                codon = seq[i:i + 3]
                if codon == "TAA":
                    break
                elif codon == "TGA":
                    break
                elif codon == "TAG":
                    break
                else:
                    frame1 += codon
    print("frame 1")
    print(frame1)
    for i in range(1, len(seq), 3):
        codon = seq[i:i + 3]
        if codon == "ATG":
            start = i
            for i in range(start, len(seq), 3):
                codon = seq[i:i + 3]
                if codon == "TAA":
                    break
                elif codon == "TGA":
                    break
                elif codon == "TAG":
                    break
                else:
                    frame2 += codon
    print("frame 2")
    print(frame2)
    for i in range(2, len(seq), 3):
        codon = seq[i:i + 3]
        if codon == "ATG":
            start = i
            for i in range(start, len(seq), 3):
                codon = seq[i:i + 3]
                if codon == "TAA":
                    break
                elif codon == "TGA":
                    break
                elif codon == "TAG":
                    break
                else:
                    frame3 += codon
    print("frame 3")
    print(frame3)
    frameList = [frame1, frame2, frame3]
    numberList = [len(frame1), len(frame2), len(frame3)]
    print(numberList)
    maxValue = max(numberList)
    maxIndex = numberList.index(maxValue)
    gene = frameList[maxIndex]
    return(gene)


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


seq = readFasta()
print(seq)
reverseComp = dnaPol(seq)
print(reverseComp)
gene = rnaPol()
print(gene)
protein = translate(gene)
print(protein)
# print(reverseComp)


# prt = read_seq("amino_acid_sequence_original.txt")
# dna = read_seq(inputfile)
# print(dna)
# p = translate(dna[20:935])
# p == prt
