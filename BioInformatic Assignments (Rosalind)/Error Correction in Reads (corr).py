# Functions


def CRI(counts, oseq):  # correct and incorrect function
    correct = []
    incorrect = []
    for s in counts:
        if counts[s] >= 2:
            correct.append(s)
        elif s in oseq:
            incorrect.append(s)
    return correct, incorrect


def CR(correct, incorrect):  # Correction Function
    corrected = []
    for seq1 in incorrect:
        for seq2 in correct:
            if HD(seq1, seq2) == 1:
                corrected.append([seq1, seq2])
    return corrected


def HD(seq1, seq2):  # Hamming Distance
    mutations = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mutations += 1
    return mutations


# Run
from Bio import SeqIO
from collections import Counter

f = open("rosalind_corr.txt", 'r')
seqs = []
oseq = []
for s in SeqIO.parse(f, "fasta"):
    oseq.append(str(s.seq))
    seqs.append(str(s.seq))
    seqs.append(str(s.seq.reverse_complement()))
f.close()
counts = Counter(seqs)
correct, incorrect = CRI(counts, oseq)
corrs = CR(correct, incorrect)
o = open("output.txt", 'w')
for i in corrs:
    print("->".join(i), file=o)

