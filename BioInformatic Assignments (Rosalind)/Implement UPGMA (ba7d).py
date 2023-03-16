
import numpy as np
from copy import deepcopy


class upgma:
    def __init__(self):
        number, dist_matrix = self.read_from_file()
        adj = self.run_upgma(dist_matrix, number)
        self.print_graph(adj)

    @staticmethod
    def read_from_file():
        file = open('rosalind_ba7d.txt', 'r')
        data = []
        for line in file:
            data.append(line.strip())
        number = int(data[0])
        dist_matrix = [[0] * number for _ in range(number)]
        for i in range(number):
            d = data[i + 1].split()
            for k in range(number):
                dist_matrix[i][k] = int(d[k])
        return number, dist_matrix

    @staticmethod
    def run_upgma(dist_matrix, number):
        dist_data = np.array(dist_matrix, dtype=float)
        np.fill_diagonal(dist_data, np.inf)
        clusters = [[i, 1] for i in range(number)]
        adj = [[] for _ in range(number)]
        age = [0. for _ in range(number)]
        if len(dist_data) <= 1:
            return adj
        while True:
            index = np.argmin(dist_data)
            i = index // len(dist_data)
            j = index % len(dist_data)
            i_new = len(adj)
            adj.append([])
            clusters_new = [i_new, clusters[i][1] + clusters[j][1]]
            adj[i_new].append(clusters[i][0])
            adj[i_new].append(clusters[j][0])
            adj[clusters[i][0]].append(i_new)
            adj[clusters[j][0]].append(i_new)
            age.append(dist_data[i, j] / 2)

            if 2 == len(dist_data):
                break

            dist_data_new = (dist_data[i, :] * clusters[i][1] + dist_data[j, :] * clusters[j][1]) / (clusters[i][1] +
                                                                                                     clusters[j][1])
            dist_data_new = np.delete(dist_data_new, [i, j], 0)
            dist_data = np.delete(dist_data, [i, j], 0)
            dist_data = np.delete(dist_data, [i, j], 1)
            dist_data = np.insert(dist_data, len(dist_data), dist_data_new, axis=0)
            dist_data_new = np.insert(dist_data_new, len(dist_data_new), np.inf, axis=0)
            dist_data = np.insert(dist_data, len(dist_data) - 1, dist_data_new, axis=1)

            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]

            clusters.append(clusters_new)

        adj_l = deepcopy(adj)
        for i, nodes in enumerate(adj):
            for j, v in enumerate(nodes):
                adj_l[i][j] = (v, abs(age[i] - age[v]))
        return adj_l

    @staticmethod
    def print_graph(adj):
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                print(str(i) + '->' + str(d) + ':' + '%0.3f' % w)


if __name__ == "__main__":
    upgma()
