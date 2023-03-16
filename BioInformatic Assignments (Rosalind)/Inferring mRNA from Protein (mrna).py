# Reading a File
with open('rosalind_mrna.txt') as file:
    dna = file.read()
    dna = dna.strip()

# Amino Acid to Codon Dictionary
ACD = {
    'A': 4,
    'C': 2,
    'D': 2,
    'E': 2,
    'F': 2,
    'G': 4,
    'H': 2,
    'I': 3,
    'K': 2,
    'L': 6,
    'M': 1,
    'N': 2,
    'P': 4,
    'Q': 2,
    'R': 6,
    'S': 6,
    'T': 4,
    'V': 4,
    'W': 1,
    'Y': 2
}

# Run
seq = list(dna)
OT = 1
for i in seq:
    OT = ACD[i] * OT

# Output
print((OT*3) % 1000000)
