from collections import defaultdict
from Bio import SeqIO
from graphviz import Digraph

from new.edge import Edge
from new.vertex import Vertex


class Graph:
    def __init__(self, path, k):
        self.reads = SeqIO.parse(path, 'fasta')
        self.k = k
        self.graph_scheme = {}
        self.edges = {}

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
            for v, out in self.graph_scheme.values():
                dot.node(v.sequence, label=f'{v.sequence} {v.coverage}')
            for e in self.edges.values():
                dot.edge(e.sequence[:-1], e.sequence[1:], label=f'{e.sequence} {self.k + len(e.coverages) - 1} {e.coverage}')
        else:
            for v, out in self.graph_scheme.values():
                dot.node(v.sequence, label=f'{v.coverage}')
            for e in self.edges.values():
                dot.edge(e.sequence[:-1], e.sequence[1:], label=f'{self.k + len(e.coverages) - 1} {e.coverage}')
        # Save pdf
        dot.render(filename, view=True)

    # def collapse(self):
    #     # If 1st vertex have 1 out degree and 2nd have 1 in degree
    #     # Add new links in graph
    #     # Delete adjacent intermediate vertex and replace links between vertex, it and next vertices, replace edges
    #     for vertex, adjacents in self.graph_scheme.values():
    #         if vertex.out_degree == 1 and self.graph_scheme[adjacents[0]][0].in_degree == 1:
    #             neighbour = adjacents[0]
    #             # Clear adjacency list
    #             first_vertex = self.graph_scheme[vertex.sequence]
    #             first_vertex[1].clear()
    #             # self.edge_set.discard(Edge(vertex, self.graph_scheme[adjacents[0]][0]))
    #             # Add new links and modify edges in set
    #             # print(neighbour)
    #             # print(self.graph_scheme[neighbour])
    #             for next_vertex in self.graph_scheme[neighbour]:
    #                 # print(first_vertex)
    #                 # print(type(first_vertex))
    #                 # print(first_vertex[1])
    #                 # print(self.graph_scheme[neighbour])
    #                 first_vertex[1].append(next_vertex[0])
    #                 for edge in self.edge_set:
    #                     if edge.source == vertex and edge.dest == self.graph_scheme[adjacents[0]][0]:
    #                         self.edge_set.discard(edge)
    #                     elif edge.source == self.graph_scheme[adjacents[0]][0] and edge.dest == next_vertex[0]:
    #                         edge.add_part(vertex)



    def fragmentate(self):
        """
        Cleave reads into kmers and add them to graph
        :return:
        """
        # Cleave all reads to kmers
        for read in self.reads:
            self.add_read(str(read.seq).upper())

    def add_read(self, read):
        # Get all sequential kmers from read
        kmers = [read[i:i + self.k] for i in range(len(read) - self.k + 1)]

        # Iterate over adjacent kmers
        for source, destination in zip(kmers, kmers[1:]):
            # Possible cases:
            # All vertices are present - increase their coverage, increase vertex degrees if there wasn't such edge and
            # add new edge
            # 1 or 2 vertices are absent - create them, increase vertex degrees and coverage, add new edge to set
            if source in self.graph_scheme and destination in self.graph_scheme:
                self.graph_scheme[source][0].increase_coverage()
                self.graph_scheme[destination][0].increase_coverage()
                if (source, destination) not in self.edges:
                    self.graph_scheme[source][0].increase_out_degree()
                    self.graph_scheme[destination][0].increase_in_degree()

                    self.graph_scheme[source][1].append(destination)

                    edge = Edge(self.graph_scheme[source][0], self.graph_scheme[destination][0])
                    self.edges[(source, destination)] = edge
                continue

            # Create absent vertex
            elif source in self.graph_scheme and destination not in self.graph_scheme:
                self.graph_scheme[destination] = (Vertex(destination), [])
            elif source not in self.graph_scheme and destination in self.graph_scheme:
                self.graph_scheme[source] = (Vertex(source), [])
            else:
                self.graph_scheme[source] = (Vertex(source), [])
                self.graph_scheme[destination] = (Vertex(destination), [])

            # Increase degrees of vertices
            self.graph_scheme[source][0].increase_out_degree()
            self.graph_scheme[destination][0].increase_in_degree()
            # Increase vertex coverage
            self.graph_scheme[source][0].increase_coverage()
            self.graph_scheme[destination][0].increase_coverage()
            # Add link in graph edge to edges set
            self.graph_scheme[source][1].append(destination)
            self.edges[(source, destination)] = Edge(self.graph_scheme[source][0], self.graph_scheme[destination][0])

    def cover_edges(self):
        for edge in self.edges.values():
            edge.ini_cover()

    def edge_coverage(self):
        for edge in self.edges.values():
            edge.compute_coverage()

    # TODO kmer all reads and create from them vertices and edges
    # TODO link adjacent kmers in reads
    # TODO run over edges, count initial coverage
    # TODO
    # TODO

a = Graph('../test4', 3)

# Get graph
# add coverage of vertices to edges
# collapse graph
# compute edge coverage
# plot graph
a.fragmentate()
print(a.graph_scheme)
a.cover_edges()

# a.collapse()
a.edge_coverage()

a.plot('fd', True)


for v in a.graph_scheme.values():
    print(v[0])

for e in a.edges.values():
    print(e)
