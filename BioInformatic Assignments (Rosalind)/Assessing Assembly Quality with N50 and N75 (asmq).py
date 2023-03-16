
def asmq(dna_list):
    total_len = sum([len(x) for x in dna_list])
    dna_list.sort(key=len)
    curr_len = 0
    n50 = total_len
    n75 = total_len
    for i in range(len(dna_list) - 1, -1, -1):
        curr_len += len(dna_list[i])
        if n50 == total_len and curr_len > total_len * 0.5:
            n50 = len(dna_list[i])
        if curr_len > total_len * 0.75:
            n75 = len(dna_list[i])
            break
    result = str(n50) + " " + str(n75)
    return result


if __name__ == "__main__":
    f = open('rosalind_asmq.txt', 'r')
    dna_list = f.read().splitlines()
    output = asmq(dna_list)
    print(output)
