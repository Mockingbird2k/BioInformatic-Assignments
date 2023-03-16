
class Node:
    def __init__(self, kmer):
        self.kmer = kmer
        self.edges = set()

    def __str__(self):
        return "kmer = {}\n edges = {}".format(self.kmer, self.edges)


def kmers_list(m, string):

    kmers = []
    for i in range(len(string) - m + 1):
        kmers.append(string[i: i + m])
    return kmers


def reverse_comp(dna_strand):

    complement_dna = ""

    dnaSequence = list(dna_strand)
    dnaSequence.reverse()

    dna_strand = ''.join(dnaSequence)
    complement_dict = {"C": "G", "G": "C", "T": "A", "A": "T"}

    for base in dna_strand:
        complement_dna += complement_dict[base]

    return complement_dna


def reverse_applier(file):
    dna_list = file.read().splitlines()
    dna_list += [reverse_comp(x) for x in dna_list]
    kmer_length = len(dna_list[0]) - 1
    return dna_list, kmer_length


if __name__ == "__main__":
    f = open('rosalind_dbru.txt', 'r')

    dna_list, kmer_length = reverse_applier(f)
    node_dict = {}

    for dna in dna_list:
        last_node = 0
        for kmer in kmers_list(kmer_length, dna):
            if kmer not in node_dict:
                node_dict[kmer] = Node(kmer)
            if last_node:
                last_node.edges.add(node_dict[kmer])
            last_node = node_dict[kmer]

    for kmer, node in sorted(node_dict.items()):
        for edge in node.edges:
            print("({0}, {1})".format(kmer, edge.kmer))
