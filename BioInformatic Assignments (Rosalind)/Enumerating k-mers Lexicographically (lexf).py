# Defining a Function

def lex(series, n, ini='', seq=[]):
    if n == 0:
        seq.append(ini)
    else:
        for c in series:
            lex(series, n - 1, ini + c, seq)
    return seq

# Run


if __name__ == "__main__":

    D = open('rosalind_lexf.txt').read()
    Data = D.split()
    series = Data[:-1]
    n = int(Data[-1])

    for i in lex(series, n):
        print(i)

