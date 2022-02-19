import os
from tracemalloc import start

# REMEMBER:
# Change dir below to choose the folder you want to read and
# write to before running this program.
os.chdir(r'C:\Users\Madison Hypes\Desktop\FOLDER\Code\ECS129Project')


def writeOutput():
    with open("SequenceOutput.txt", "w") as f:
        f.write(
            '\n.-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-.')
        f.write('\nWelcome to our ECS129 project for 2022!\n')
        f.write('\n')
        f.write('Project made by Madison Hypes, Jason Hu')
        f.write(', and David Galkowski.\n')
        f.write('Created on February 19th, 2022\n')
        f.write('\n')
        f.write('Our project will read a dna sequence in')
        f.write(' fasta format, and return\nthe protein made')
        f.write(' from the longest open reading frame.\n')
        f.write('\nRESULTS SHOWN BELOW\n-------\n')
        f.write('Input Sequence:\n' + seq + '\n')
        f.write('\nReverse Complement Sequence:\n')
        f.write(reverseComp)
        f.write('\n\nLargest ORF in Input Sequence:\n')
        f.write(globalList[0])
        f.write('\n\nLargest ORF in Reverse Complement ')
        f.write('Sequence:\n')
        f.write(globalList[1] + '\n\n')
        f.write('LARGEST ORF between Input Sequence and ')
        f.write('Reverse Complement Sequence -\n')
        if len(globalList[0]) > len(globalList[1]):
            f.write('ORF from Input Sequence is Larger')
            f.write(' and is written below:\n')
            f.write(globalList[0])
        elif len(globalList[0]) == len(globalList[1]):
            f.write('ORF from Inputs Sequence is equal length')
            f.write(' to ORF from Reverse Complement:\n')
            f.write(globalList[0])
        else:
            f.write('ORF from Reverse Complement is Larger')
            f.write(' and is written below:\n')
            f.write(globalList[1])
        f.write('\n\nmRNA Sequence of ORF:\n')
        f.write(mRNASeq)
        f.write('\n\nProtein sequence of ORF:\n')
        f.write(protein)
        f.write('\n\nPOSSIBLE EXTENSIONS\n')
        f.write('-------\n')
        f.write('1.) Create pipeline from ORF to homology')
        f.write(' modeling software like SWISSMODEL.\n')
        f.write('2.) Create pipeline searching AlphaFold with ORF.\n')
        f.write('3.) Connect a thermostabilty engineering program ')
        f.write('using BLOSOM 62 and swaping for thermostable\n')
        f.write('amino acid substitutes.\n\n')
        f.write('Thank you for visiting our program. Have a')
        f.write(' great day!')

        f.write('\n`-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        f.write('-=-=-=-=-=-=-=-=-=-=-=-=-=')
        f.write('-=-=-=-=-=-=-=-=-´')
        f.close()


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


def rnaPol(seq):
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
            break
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
            break
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
            break
    frameList = [frame1, frame2, frame3]
    numberList = [len(frame1), len(frame2), len(frame3)]
    maxValue = max(numberList)
    maxIndex = numberList.index(maxValue)
    gene = frameList[maxIndex]
    globalList.append(gene)
    return(globalList)


def longestORF(globalList):
    numberList = []
    numberList = [len(globalList[0]), len(globalList[1])]
    maxValue = max(numberList)
    maxIndex = numberList.index(maxValue)
    ORF = globalList[maxIndex]
    return(ORF)


def mRNA(ORF):
    table = {"A": "A", "T": "U", "G": "G", "C": "C"}
    mRNASeq = ""
    for i in range(0, len(ORF), 1):
        bp = ORF[i]
        mRNASeq += table[bp]
    return mRNASeq


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


globalList = []
seq = readFasta()
reverseComp = dnaPol(seq)
rnaPol(seq)
rnaPol(reverseComp)
ORF = longestORF(globalList)
mRNASeq = mRNA(ORF)
protein = translate(ORF)
writeOutput()

# Thank you for using our program.
# hypesmadison@gmail.com
