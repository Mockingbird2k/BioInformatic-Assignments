def genome_assembly(edgelist):
    elem = edgelist[0]
    superstring = elem[0][-1]
    i = 0
    next = 0
    while i < len(edgelist) - 1:
        superstring += elem[1][-1]
        for j in range(len(edgelist)):
            if elem[1] == edgelist[j][0]:
                next = j
                break
        elem = edgelist[next]
        i += 1
    return superstring


def construct_graph(inputList):
    edgelist = []
    length = len(inputList[0])
    inputSet = set(inputList)
    for item in inputSet:
        left = item[0:length - 1]
        right = item[1:length]
        edgelist.append([left, right])
    return edgelist


if __name__ == "__main__":
    f = open('rosalind_pcov.txt', 'r')
    dna_list = f.read().splitlines()
    new_set = construct_graph(dna_list)
    output = genome_assembly(new_set)
    print(output)
