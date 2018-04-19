from collections import defaultdict
from Bio import SeqIO
from vertex import Vertex
from edge import Edge
from graphviz import Digraph


class GraphDB:
    def __init__(self, path, k):
        """
        GraphDB initializer
        :param path: str - path to fastq file
        :param k: int - kmer size
        """
        self.path = path
        self.k = k
        # Graph representation with Vertex as keys and set of Edges as values
        self.graph = defaultdict(set)

    def add_read(self, read):
        """
        Fragment read to kmers and add them to graph
        :param read:
        :return:
        """
        # Get all sequential kmers from read
        kmers = [read[i:i + self.k] for i in range(len(read) - self.k + 1)]

        # For all kmers and their shift by 1 nucleotide (next kmer):
        # Create Vertex, add new Edge to set in graph dictionary
        # Increase coverage of Vertex corresponding to 1st kmer
        # I`ve made storing of vertices coverage in Edge because I don`t know how to obtain appropriate
        # vertex after dictionary look-up beatifully.
        # Note the last kmer didn`t get Vertex object and coverage, it will obtain it from other reads
        # if they have such kmer
        for source, destination in zip(kmers, kmers[1:]):
            e = Edge(source[0] + destination)
            self.graph[Vertex(source)].add(e)
            # next(iter(self.graph[Vertex(source)])).increase_source_coverage()
            for edge in self.graph[Vertex(source)]:
                if edge == e:
                    edge.increase_source_coverage()

    def compute_coverage(self):
        """
        Computes vertex coverage
        :return:
        """
        for vertex, edges in self.graph.items():
            vertex.coverage = str(sum(edge.source_coverage for edge in edges))
            # TODO remove source_coverage and assign it to coverage
            for edge in edges:
                edge.coverage = edge.source_coverage

    # def compute_edge_coverage(self):
    #     for edges in self.graph.values():
    #         for edge in edges:
    #             edge.coverage = self.graph[Vertex(edge.sequence[:-1])]

    def fragmentate(self):
        """
        Cleave reads into kmers and add them to graph
        :return:
        """
        # Read data
        reads = SeqIO.parse(self.path, 'fasta')
        # Cleave all reads to kmers
        for read in reads:
            self.add_read(str(read.seq).upper())
            # self.add_read(str(read.reverse_complement().seq).upper())

    def plot(self, filename, include_seq=False):
        """
        Plot assembly graph
        :param filename: str - path to output
        :param include_seq: boolean - whether to include sequences in vertices and edges
        :return:
        """
        dot = Digraph(comment='Assembly')

        # Iterate over vertices and edges, add them to graphviz graph.
        # There will be coverages of vertices, edges and their length and corresponding sequences if include_seq
        if include_seq:
            for vertex in self.graph:
                dot.node(vertex.sequence, label=f'{vertex.sequence} {vertex.coverage}')
            for edges in self.graph.values():
                for edge in edges:
                    dot.edge(edge.sequence[:-1], edge.sequence[1:], label=f'{edge.sequence} {self.k + edge.size - 1} {edge.coverage}')
        else:
            for vertex in self.graph:
                dot.node(vertex.sequence, label=f'{vertex.coverage}')
            for edges in self.graph.values():
                for edge in edges:
                    dot.edge(edge.sequence[:-1], edge.sequence[1:], label=f'{self.k + edge.size - 1} {edge.coverage}')
        # Save pdf
        dot.render(filename, view=True)


if __name__ == '__main__':
    a = GraphDB('test4', 3)
    a.fragmentate()
    a.compute_coverage()

    for vert, edges in a.graph.items():
        print(f'{vert}, {list(map(str, edges))}')

    a.plot('c', True)
