# Load libraries
import numpy as np

# Main
Fasta = np.genfromtxt("rosalind_gc.txt", dtype=str)  # Data Input
m, c, g, leng, maxx, k, maxi = 0, 0, 0, 0, 0, 0, 0  # Initial Values
for i in range(len(Fasta)):
    # Sequence Calculator
    if Fasta[i][0] != '>':
        leng += len(Fasta[i])
        c += np.char.count(Fasta[i], 'C')
        g += np.char.count(Fasta[i], 'G')

    # Max Calculation
    if i > 0 and Fasta[i][0] == '>' or i == len(Fasta)-1:
        if maxx < ((c + g) / leng) * 100:
            maxx = ((c + g) / leng) * 100
            maxi = k
        c, g, leng = 0, 0, 0

    # Title Index
    if Fasta[i][0] == '>':
        k = i
# Results
print(Fasta[maxi][1:])
print(maxx)
