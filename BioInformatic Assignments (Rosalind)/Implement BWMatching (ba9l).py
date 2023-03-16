
def preprocess_bwt(bwt, alphabet):

    l = len(bwt)
    counts = dict()
    starts = dict()
    for char in alphabet:
        counts[char] = [0] * (l + 1)
    for i in range(l):
        currChar = bwt[i]
        for char, count in counts.items():
            counts[char][i + 1] = counts[char][i]
        counts[currChar][i + 1] += 1
    currIndex = 0
    for char in sorted(alphabet):
        starts[char] = currIndex
        currIndex += counts[char][l]
    return starts, counts


def count_occurrences(pattern, bwt, starts, counts):

    top = 0
    bottom = len(bwt) - 1
    currIndex = len(pattern) - 1
    while top <= bottom:
        if currIndex >= 0:
            symbol = pattern[currIndex]
            currIndex -= 1
            if counts[symbol][bottom + 1] - counts[symbol][top] > 0:
                top = starts[symbol] + counts[symbol][top]
                bottom = starts[symbol] + counts[symbol][bottom + 1] - 1
            else:
                return 0
        else:
            return bottom - top + 1


if __name__ == '__main__':

    f = open('rosalind_ba9l.txt', 'r')
    dna_list = f.read().strip().split()
    bwt = dna_list[0]
    patterns = dna_list[1:]
    alphabet = ['$', 'A', 'C', 'G', 'T']
    starts, occ_counts_before = preprocess_bwt(bwt, alphabet)
    occurrence_counts = []
    for pattern in patterns:
        occurrence_counts.append(count_occurrences(pattern, bwt, starts, occ_counts_before))
    print(' '.join(map(str, occurrence_counts)))
