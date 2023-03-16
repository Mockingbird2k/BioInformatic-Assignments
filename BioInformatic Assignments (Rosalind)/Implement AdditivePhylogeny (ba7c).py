
import queue


class additive_phylogeny:

    def __init__(self):
        n, dist_matrix = self.read_from_file()
        adj = self.reconstructPhylogeny(dist_matrix, n)
        self.printGraph(adj)

    def read_from_file(self):
        file = open('rosalind_ba7c.txt', 'r')
        data = []
        for line in file:
            data.append(line.strip())
        n = int(data[0])
        dist_matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            d = data[i + 1].split()
            for k in range(n):
                dist_matrix[i][k] = int(d[k])
        return n, dist_matrix

    def calculate_limb_lenth(self, dist_matrix, n, j):
        limb_length = float('inf')
        if j > 0:
            i = j - 1
        else:
            i = j + 1
        for k in range(n):
            if i != k and k != j:
                curr_length = (dist_matrix[i][j] + dist_matrix[j][k] - dist_matrix[i][k]) // 2
                if curr_length < limb_length:
                    limb_length = curr_length
                    curr_index = (i, k)
        return limb_length, curr_index[0], curr_index[1]

    def reconstructPhylogeny(self, data, n):

        def addNode(adj, j, limb_length, i, k, x):
            l = len(adj)
            dist = [float('inf')] * l
            parent = [-1] * l
            q = queue.Queue()
            dist[i] = 0
            q.put(i)
            while not q.empty():
                curr_node = q.get()
                for node, weight in adj[curr_node].items():
                    if float('inf') == dist[node]:
                        dist[node] = dist[curr_node] + weight
                        parent[node] = curr_node
                        q.put(node)
                        if node == k:
                            prevNode = node
                            while x < dist[prevNode]:
                                curr_node = prevNode
                                prevNode = parent[curr_node]
                            if x == dist[prevNode]:
                                adj[prevNode][j] = limb_length
                                adj[j][prevNode] = limb_length
                            else:
                                adj.append(dict())
                                newNode = len(adj) - 1
                                adj[j][newNode] = limb_length
                                adj[newNode][j] = limb_length
                                del adj[prevNode][curr_node]
                                del adj[curr_node][prevNode]
                                adj[prevNode][newNode] = x - dist[prevNode]
                                adj[newNode][prevNode] = x - dist[prevNode]
                                adj[curr_node][newNode] = dist[curr_node] - x
                                adj[newNode][curr_node] = dist[curr_node] - x
                            return

        adj = [dict() for _ in range(n)]
        adj[0][1] = data[0][1]
        adj[1][0] = data[1][0]
        for j in range(2, n):
            limb_length, i, k = self.calculate_limb_lenth(data, j + 1, j)
            x = data[i][j] - limb_length
            addNode(adj, j, limb_length, i, k, x)
        return adj

    def printGraph(self, adj):
        for i, dicts in enumerate(adj):
            for d, w in dicts.items():
                print(str(i) + '->' + str(d) + ':' + str(w))


if __name__ == "__main__":
    additive_phylogeny()
