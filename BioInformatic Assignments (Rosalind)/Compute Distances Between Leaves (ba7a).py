import queue


class distance_between_leaves:

    def __init__(self):

        n, adj_dict = self.read_file()
        dist_matrix = self.calculate_distance_matrix(n, adj_dict)
        self.print_distance_matrix(dist_matrix)

    def read_file(self):

        file = open('rosalind_ba7a.txt', 'r')
        data = []
        for line in file:
            data.append(line.strip())
        number = int(data[0])
        adj_dict = dict()
        for d in data[1:]:
            d = d.split('->')
            d1 = d[1].split(':')
            if not int(d[0]) in adj_dict:
                adj_dict[int(d[0])] = []
            adj_dict[int(d[0])].append((int(d1[0]), int(d1[1])))

        return number, adj_dict

    def calculate_distance_matrix(self, number, adj_dict):

        dist_matrix = [[0] * number for _ in range(number)]
        for i in range(number):
            dist = dict()
            q = queue.Queue()
            dist[i] = 0
            q.put(i)
            while not q.empty():
                curr_node = q.get()
                for node, weight in adj_dict[curr_node]:
                    if not node in dist:
                        dist[node] = dist[curr_node] + weight
                        if node < number:
                            dist_matrix[i][node] = dist[node]
                        q.put(node)

        return dist_matrix

    def print_distance_matrix(self, dist_matrix):

        for j in dist_matrix:
            print(' '.join([str(i) for i in j]))


if __name__ == "__main__":
    distance_between_leaves()
