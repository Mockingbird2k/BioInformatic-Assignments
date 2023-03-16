
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
with open('rosalind_pdst.txt') as file:
    dna_list = fasta(file)

dna_list = [x.content for x in dna_list]
n = len(dna_list)
s_len = len(dna_list[0])
assert all([len(x) == s_len for x in dna_list])
D = [[0.0 for j in range(n)] for i in range(n)]
for i in range(n - 1):
    for j in range(i + 1, n):
        D[i][j] = D[j][i] = sum([dna_list[i][k] != dna_list[j][k] for k in range(s_len)]) / s_len

for row in D:
    print(' '.join(['{:.5f}'.format(x) for x in row]))
