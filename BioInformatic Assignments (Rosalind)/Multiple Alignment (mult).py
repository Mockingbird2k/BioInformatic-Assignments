# Class

class DNA:
    def __init__(self, name='', info='', content=None):
        self.name = name
        self.info = info
        self.content = content

    def __repr__(self):
        return "name: {}\ninfo: {}\ncontent: {}".format(self.name, self.info, self.content)

# Functions


def anD(dna_list, line):  # anD = add new DNA
    assert line[0] == '>'
    first_space_idx = line.find(' ')
    if first_space_idx != -1:
        dna_name = line[1:first_space_idx]
        dna_info = line[first_space_idx:].strip()
    else:
        dna_name = line[1:]
        dna_info = ''
    dna_list.append(DNA(name=dna_name, info=dna_info, content=[]))


def alD(cur_DNA, line):  # alD = add line DNA
    for x in line:
        if x in VC:
            cur_DNA.content.append(x)
        elif x == ' ':
            continue
        else:
            raise Exception()


def fasta(file):

    state = 0
    dna_list = []
    for line in file:
        line = line.strip()
        if state == 0:
            if line[0] == '>':
                anD(dna_list, line)
                state = 1
            elif line == '':
                continue
            else:
                raise Exception()
        elif state == 1:
            alD(dna_list[-1], line)
            state = 2
        elif state == 2:
            if line[0] == '>':
                anD(dna_list, line)
                state = 1
            else:
                alD(dna_list[-1], line)
        else:
            raise Exception()
    file.seek(0)
    return dna_list


# Run
VC = {'A', 'T', 'C', 'G'}
with open('rosalind_mult.txt') as file:
    dna_list = fasta(file)

a, b, c, d = [x.content for x in dna_list]


def get_best_score(D, i, j, k, l, ax, bx, cx, dx, match_cost, mismatch_cost, sign):

    perms = [
        [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, 0, 0], [0, 1, 0, 1],
        [0, 1, 1, 0], [0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, 0],
        [1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 1, 1, 1]
    ]
    best_score = -99999
    best_perm = None
    for perm in perms:
        # Check if perm is applicable.
        if not all([not perm[0] or i > 0, not perm[1] or j > 0, not perm[2] or k > 0, not perm[3] or l > 0]):
            continue
        ax_p = ax if perm[0] else '-'
        bx_p = bx if perm[1] else '-'
        cx_p = cx if perm[2] else '-'
        dx_p = dx if perm[3] else '-'
        tmp = [ax_p, bx_p, cx_p, dx_p]
        score = D[i - perm[0]][j - perm[1]][k - perm[2]][l - perm[3]]  # Include corresponding previous D value.

        for p0 in range(len(tmp) - 1):
            for p1 in range(p0 + 1, len(tmp)):
                score += sign * match_cost if tmp[p0] == tmp[p1] else sign * mismatch_cost

        if score > best_score:
            best_score = score
            best_perm = perm
    return best_score, best_perm


match_cost = 0
mismatch_cost = 1
sign = -1
D = [[[[0 for l in range(len(d) + 1)] for k in range(len(c) + 1)] for j in range(len(b) + 1)] for i in
     range(len(a) + 1)]
backtrace = [[[[None for l in range(len(d) + 1)] for k in range(len(c) + 1)] for j in range(len(b) + 1)] for i in
             range(len(a) + 1)]
for i in range(1, len(a) + 1):
    D[i][0][0][0] = i * sign * 3 * mismatch_cost
    backtrace[i][0][0][0] = [1, 0, 0, 0]
for j in range(1, len(b) + 1):
    D[0][j][0][0] = j * sign * 3 * mismatch_cost
    backtrace[0][j][0][0] = [0, 1, 0, 0]
for k in range(1, len(c) + 1):
    D[0][0][k][0] = k * sign * 3 * mismatch_cost
    backtrace[0][0][k][0] = [0, 0, 1, 0]
for l in range(1, len(d) + 1):
    D[0][0][0][l] = l * sign * 3 * mismatch_cost
    backtrace[0][0][0][l] = [0, 0, 0, 1]
for i in range(len(a) + 1):
    for j in range(len(b) + 1):
        for k in range(len(c) + 1):
            for l in range(len(d) + 1):
                # If at one of the base case lines, skip.
                if i + j + k + l <= 1:
                    continue
                best_score, best_perm = get_best_score(D, i, j, k, l, a[i - 1], b[j - 1], c[k - 1], d[l - 1],
                                                       match_cost, mismatch_cost, sign)
                D[i][j][k][l] = best_score
                backtrace[i][j][k][l] = best_perm

i, j, k, l = len(a), len(b), len(c), len(d)
a_align = []
b_align = []
c_align = []
d_align = []
while i > 0 or j > 0 or k > 0 or l > 0:
    perm = backtrace[i][j][k][l]
    a_align.append(a[i - 1] if perm[0] else '-')
    b_align.append(b[j - 1] if perm[1] else '-')
    c_align.append(c[k - 1] if perm[2] else '-')
    d_align.append(d[l - 1] if perm[3] else '-')
    i, j, k, l = i - perm[0], j - perm[1], k - perm[2], l - perm[3]

print(D[-1][-1][-1][-1])
print(''.join(a_align[::-1]))
print(''.join(b_align[::-1]))
print(''.join(c_align[::-1]))
print(''.join(d_align[::-1]))

# Verify alignment score.
alignment_score = 0
tmp = [a_align, b_align, c_align, d_align]
for i in range(len(a_align)):
    for j in range(len(tmp) - 1):
        for k in range(j + 1, len(tmp)):
            if tmp[j][i] != tmp[k][i]:
                alignment_score -= 1
error_msg = "D value -> {} != {} <- score of result alignment".format(D[-1][-1][-1][-1], alignment_score)
assert D[-1][-1][-1][-1] == alignment_score, error_msg
