# Load libraries
import numpy as np

# Main

with open("rosalind_revc.txt", 'r') as d:
    data = d.readline()
s = np.zeros(len(data)-1, dtype=str)
for i in range(len(s)):
    # Sequence Calculator
    if data[i] == 'A':
        s[len(s)-i-1] = 'T'
    elif data[i] == 'T':
        s[len(s)-i-1] = 'A'
    elif data[i] == 'C':
        s[len(s)-i-1] = 'G'
    elif data[i] == 'G':
        s[len(s)-i-1] = 'C'
for i in range(len(s)):
    if i < len(s)-1:
        print(s[i], end='')
    else:
        print(s[i])

