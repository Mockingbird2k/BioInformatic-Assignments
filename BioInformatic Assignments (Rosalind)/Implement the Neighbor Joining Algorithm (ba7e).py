
import numpy as np


class imp_neighbor_joining:
    def __init__(self):
        number, dist_matrix = self.read_from_file()
        adj = self.run_neighbor_joining(dist_matrix, number)
        self.print_graph(adj)

    @staticmethod
    def read_from_file():
        file = open('rosalind_ba7e.txt', 'r')
        data = []
        for line in file:
            data.append(line.strip())
        number = int(data[0])
        dist_matrix = [[0] * number for _ in range(number)]
        for i in range(number):
            d = data[i + 1].split()
            for j in range(number):
                dist_matrix[i][j] = int(d[j])
        return number, dist_matrix

    @staticmethod
    def run_neighbor_joining(dist_matrix, number):
        dist_data = np.array(dist_matrix, dtype=float)
        clusters = [i for i in range(number)]
        adj = [[] for _ in range(number)]
        if len(dist_data) <= 1:
            return adj
        while 0 == 0:
            if 2 == number:
                adj[len(adj) - 1].append((len(adj) - 2, dist_data[0][1]))
                adj[len(adj) - 2].append((len(adj) - 1, dist_data[0][1]))
                break
            total_dist = np.sum(dist_data, axis=0)
            d1 = (number - 2) * dist_data - total_dist - total_dist.reshape((number, 1))
            np.fill_diagonal(d1, 0.)
            index = np.argmin(d1)
            ind_number_fratio = index // number
            ind_number_mod = index % number
            delta = (total_dist[ind_number_fratio] - total_dist[ind_number_mod]) / (number - 2)
            l_ind_number_fratio = (dist_data[ind_number_fratio, ind_number_mod] + delta) / 2
            l_ind_number_mod = (dist_data[ind_number_fratio, ind_number_mod] - delta) / 2
            d_new = (dist_data[ind_number_fratio, :] + dist_data[ind_number_mod, :] - dist_data[ind_number_fratio,
                     ind_number_mod]) / 2
            dist_data = np.insert(dist_data, number, d_new, axis=0)
            d_new = np.insert(d_new, number, 0., axis=0)
            dist_data = np.insert(dist_data, number, d_new, axis=1)
            dist_data = np.delete(dist_data, [ind_number_fratio, ind_number_mod], 0)
            dist_data = np.delete(dist_data, [ind_number_fratio, ind_number_mod], 1)
            length_adj = len(adj)
            adj.append([])
            adj[length_adj].append((clusters[ind_number_fratio], l_ind_number_fratio))
            adj[clusters[ind_number_fratio]].append((length_adj, l_ind_number_fratio))
            adj[length_adj].append((clusters[ind_number_mod], l_ind_number_mod))
            adj[clusters[ind_number_mod]].append((length_adj, l_ind_number_mod))
            if ind_number_fratio < ind_number_mod:
                del clusters[ind_number_mod]
                del clusters[ind_number_fratio]
            else:
                del clusters[ind_number_fratio]
                del clusters[ind_number_mod]
            clusters.append(length_adj)
            number -= 1

        return adj

    @staticmethod
    def print_graph(adj):
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                print(str(i) + '->' + str(d) + ':' + '%0.3f' % w)


if __name__ == "__main__":
    imp_neighbor_joining()
