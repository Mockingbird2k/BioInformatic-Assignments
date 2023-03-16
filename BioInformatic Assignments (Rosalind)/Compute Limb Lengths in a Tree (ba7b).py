
class limb_length:
    def __init__(self):
        number, number_1, dist_matrix = self.read_from_file()
        length_limb = self.calculate_limb_lenth(number, number_1, dist_matrix)
        print(length_limb)

    def read_from_file(self):
        file = open('rosalind_ba7b.txt', 'r')
        data = []
        for line in file:
            data.append(line.strip())
        number = int(data[0])
        number_1 = int(data[1])
        dist_matrix = [[0] * number for _ in range(number)]
        for i in range(number):
            d = data[i + 2].split()
            for j in range(number):
                dist_matrix[i][j] = int(d[j])
        return number, number_1, dist_matrix

    def calculate_limb_lenth(self, n, number_1, dist_matrix):
        length_limb = float('inf')
        if number_1 > 0:
            i = number_1 - 1
        else:
            i = number_1 + 1
        for j in range(n):
            if i != j and j != number_1:
                curr_length = (dist_matrix[i][number_1] + dist_matrix[number_1][j] - dist_matrix[i][j]) // 2
                if curr_length < length_limb:
                    length_limb = curr_length
        return length_limb


if __name__ == "__main__":
    limb_length()
