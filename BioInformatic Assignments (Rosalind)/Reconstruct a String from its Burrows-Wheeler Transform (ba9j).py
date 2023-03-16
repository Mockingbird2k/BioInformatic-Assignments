
def inverse_bwt(bwt):

    len_bwt = len(bwt)
    out = [''] * len_bwt
    count = dict()
    s_cut = [0] * len_bwt
    strt = dict()
    s_cut, count = reverse_applier(bwt, s_cut, count)
    curr_index = 0
    first_col = []
    for char, curr_count in sorted(count.items(), key=lambda x: x[0]):
        first_col += [char] * curr_count
        strt[char] = curr_index
        curr_index += curr_count
    curr_index = 0
    out = re_reverse(len_bwt, out, first_col, curr_index, strt, s_cut, bwt)

    return ''.join(out)


def reverse_applier(bwt, shortcut, count):

    l = len(bwt)
    for i in range(l):
        last_char = bwt[i]
        curr_count = count.get(last_char, 0)
        shortcut[i] = curr_count
        count[last_char] = curr_count + 1

    return shortcut, count


def re_reverse(len_bwt, out, first_col, curr_index, strt, s_cut, bwt):

    for i in range(len_bwt):
        out[len_bwt - i - 1] = first_col[curr_index]
        curr_index = strt[bwt[curr_index]] + s_cut[curr_index]

    return out


if __name__ == '__main__':
    f = open('rosalind_ba9j.txt', 'r')
    data = f.read().strip()
    output = inverse_bwt(data)
    print(output)
