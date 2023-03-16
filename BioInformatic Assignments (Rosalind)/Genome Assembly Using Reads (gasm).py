
def reverse_comp(dna_strand):

    complement_dna = ""

    dnaSequence = list(dna_strand)
    dnaSequence.reverse()

    dna_strand = ''.join(dnaSequence)
    complement_dict = {"C": "G", "G": "C", "T": "A", "A": "T"}

    for base in dna_strand:
        complement_dna += complement_dict[base]

    return complement_dna


if __name__ == '__main__':
    dna_strand = [x.strip() for x in open('rosalind_gasm.txt').readlines()]

    for k in reversed(range(2, len(dna_strand[0]) + 1)):
        mers = set()

        for s in dna_strand:
            for i in range(len(s) - k + 1):
                mers.add(s[i:i + k])
                mers.add(reverse_comp(s[i:i + k]))

        ans = None

        stack = [([], set())]

        while len(stack) > 0 and ans is None:
            path, vs = stack.pop()

            if len(path) == 0:
                for mer in mers:
                    stack.append(([mer], {mer}))
            else:
                mer = path[-1]

                for a in "A, C, T, G":
                    nmer = mer[1:] + a

                    if nmer in mers and nmer != reverse_comp(mer):
                        if nmer == path[0]:
                            ans = list(path)
                            break
                        elif nmer not in vs:
                            stack.append((path + [nmer], vs.union({nmer})))

        if ans is not None:
            output = ans[0]

            for i in range(1, len(ans)):
                output += ans[i][-1]

            output = output[:-(k - 1)]
            doutput = output + output
            success = True

            for dna in dna_strand:
                if doutput.find(dna) == -1 and doutput.find(reverse_comp(dna)) == -1:
                    success = False

            if success:
                print(output)
                break
