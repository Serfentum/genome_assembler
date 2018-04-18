from collections import defaultdict
from Bio import SeqIO
from vertex import Vertex
from edge import Edge


class GraphDB:
    # Make reverse complement of reads


    def __init__(self, path, k):
        # Dict with keys=vertex and values=[edges]
        self.path = path
        self.graph = defaultdict(list)
        self.k = k

    def add_read(self, read):
        """
        Fragment read to kmers and add them to graph
        :param read:
        :return:
        """
        kmers = [read[i:i + self.k] for i in range(len(read) - self.k + 1)]

        for source, destination in zip(kmers, kmers[1:]):
            # self.vertices[Vertex(destination, source)] = Vertex(latter)
            self.graph[Vertex(source)].append(Vertex(destination))

    def fragmentate(self):
        reads = SeqIO.parse(self.path, 'fasta')

        for read in reads:
            self.add_read(read)

    def compute_edge_coverage(self):
        # Go through graph and compute coverage of all edges given vertex coverage

    def add_edge(self, source, destination):
        edge = Edge(source.sequence + destination.sequence[-1])
        source.out_edges[destination].append(edge)
        destination.in_edges[source].append(edge)


