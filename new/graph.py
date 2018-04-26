from collections import defaultdict
from Bio import SeqIO
from graphviz import Digraph
import logging
import time

from new.edge import Edge
from new.vertex import Vertex


logging.basicConfig(filename='graph.log',
                    filemode='w',
                    format='%(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)



class Graph:
    def __init__(self, path, k):
        self.reads = SeqIO.parse(path, 'fasta')
        self.k = k
        self.graph_scheme = {}
        self.edges = {}
        self.collapsed_graph = {}
        self.collapsed_edges = {}


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
            print('EDGES')
            for e in self.edges.values():
                print(e)
            for vs, (v, out) in self.graph_scheme.items():
                dot.node(vs, label=f'{vs} {v.coverage}')
                for dv in out:
                    dot.edge(vs, dv, label=f'{self.edges[(vs, dv)].sequence} {self.k + len(self.edges[(vs, dv)].coverages) - 1} {self.edges[(vs, dv)].coverage}')
        else:
            for vs, (v, out) in self.graph_scheme.items():
                dot.node(vs, label=f'{v.coverage}')
                for dv in out:
                    dot.edge(vs, dv, label=f'{self.k + len(self.edges[(vs, dv)].coverages) - 1} {self.edges[(vs, dv)].coverage}')
        # Save pdf
        dot.render(filename, view=True)

    def collapse(self, pause=1):

        obsolete = set()
        # i = 1
        # Iterate over vertices in current graph
        logging.debug('Start collapsing:')
        for source_seq, (source, adj_vertex) in self.graph_scheme.items():
            logger.debug('%s\n\t%s is a source', obsolete, source.sequence)
            # If source vertex was marked as deleted, skip it
            if source_seq in obsolete:
                logger.debug('\t%s is in obsolete - skip', source.sequence)
                continue
            # If source vertex is a leaf or internal vertex doesn't have 1 in and 1 out edges - add them to new graph
            # and continue
            if not adj_vertex:
                logger.debug('\t%s has no adjacent vertex - add_unchanged', source.sequence)
                self.add_unchanged(source)

            for inter in adj_vertex:  # inter = self.graph_scheme[adj_vertex[0]][0] # for inter in self.graph_scheme[adj_vertex[0]]
                logger.debug('\t\tIterate over intermediate vertices')
                inter = self.graph_scheme[inter][0]
                if inter.sequence in obsolete:
                    logger.debug('\t\tInter %s is in obsolete - skip', inter.sequence)
                    continue

                if inter.in_degree != 1 or inter.out_degree != 1:
                    logger.debug('\t\tInter %s degrees are not singular - add_unchanged(%s, %s)', inter.sequence, source.sequence, inter.sequence)
                    self.add_unchanged(source, inter)
                    continue

                # List of all 3rd layer vertices, adjacent to inter vertex
                dests = [self.graph_scheme[dest][0] for dest in self.graph_scheme[inter.sequence][1]]
                # Add new links and edges in collapsed graph
                for dest in dests:
                    logger.debug('\t\t\tSource, Inter, Dest - %s %s %s', source.sequence, inter.sequence, dest.sequence)
                    if dest.sequence in obsolete:
                        logger.debug('\t\t\t\t%s is obsolete', dest.sequence)

                        if source not in self.collapsed_graph:
                            logger.debug('\t\t\t\tAdd source %s', source.sequence)
                            self.collapsed_graph[source.sequence] = source, []
                        self.collapsed_graph[source.sequence][1].append(inter.sequence)
                        logger.debug('\t\t\t\tAdd link from source %s to inter %s', source.sequence, inter.sequence)
                        self.collapsed_edges[(source.sequence, inter.sequence)] = self.edges[(source.sequence, inter.sequence)]
                        logger.debug('\t\t\t\tAdd unchanged edge between %s %s', source.sequence, inter.sequence)
                        # print(f'\tAdd existed edge between {source.sequence} and {inter.sequence}')
                        continue
                    logger.debug('\t\t\t\tExtend edge with source, inter, dest %s %s %s', source.sequence, inter.sequence, dest.sequence)
                    self.extend_one_edge(source, inter, dest)
                    # print(f'\tAdd edge between {source.sequence} and {dest.sequence}')
                logger.debug('\t\tAdd %s to obsolete', inter.sequence)
                obsolete.add(inter.sequence)
                self.slideshow(pause)
                # print(f'Discard {inter.sequence}\n')

    def slideshow(self, pause):
        dot = Digraph(comment='Assembly')
        for vs, (v, out) in self.collapsed_graph.items():
            dot.node(vs, label=f'{vs} {v.coverage}')
            for dv in out:
                dot.edge(vs, dv, label=f'{self.collapsed_edges[(vs, dv)].sequence} {self.k + len(self.collapsed_edges[(vs, dv)].coverages) - 1} {self.collapsed_edges[(vs, dv)].coverage}')
        dot.render(view=True)
        time.sleep(pause)

    def extend_one_edge(self, source, inter, dest):
        # Add source vertex to collapsed graph
        if source.sequence not in self.collapsed_graph:
            self.collapsed_graph[source.sequence] = source, []
        # If inter vertex was in collapsed graph delete it
        if inter.sequence in self.collapsed_graph:
            del self.collapsed_graph[inter.sequence]
            # del self.collapsed_edges[(inter.sequence, dest.sequence)]
            logger.debug('%s', self.collapsed_graph.keys())
        # Add link between 1st and 3rd vertices
        self.collapsed_graph[source.sequence][1].append(dest.sequence)
        # Add edge to collapsed edges
        self.modify_edge(source, inter, dest)

    def modify_edge(self, source, inter, dest):
        """
        Assign edge between source and dest vertices in collapsed graph, update edge's sequence, coverages and dest
        :param source: Vertex - 1st layer vertex
        :param inter: Vertex - 2nd layer vertex, next to source
        :param dest: Vertex - 3rd layer vertex, next to inter
        :return:
        """
        self.collapsed_edges[(source.sequence, dest.sequence)] = self.edges[(source.sequence, inter.sequence)]
        self.collapsed_edges[(source.sequence, dest.sequence)].extend(self.edges[(inter.sequence, dest.sequence)])

    def add_unchanged(self, source, inter=None):
        """
        Add source vertex or source vertex with inter and their edge into new graph intact
        :param source: Vertex - first vertex
        :param inter: Vertex - second vertex
        :return:
        """
        if source.sequence not in self.collapsed_graph:
            self.collapsed_graph[source.sequence] = source, []
        if inter and inter.sequence not in self.collapsed_graph:
            self.collapsed_graph[inter.sequence] = inter, []

        if inter:
            self.collapsed_graph[source.sequence][1].append(inter.sequence)
            # self.collapsed_graph[inter.sequence] = self.graph_scheme[inter.sequence]
            self.collapsed_edges[(source.sequence, inter.sequence)] = self.edges[(source.sequence, inter.sequence)]



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

        # Iterate over adjacent kmers, add them into graph dict as vertices, add edges between them in the edge graph
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
        # Write information about node coverage into edges
        for edge in self.edges.values():
            edge.ini_cover()

    def edge_coverage(self):
        # Compute coverage of edges
        for edge in self.edges.values():
            edge.compute_coverage()


a = Graph('../test4', 3)
# a = Graph('../wloop', 3)

# Get graph
# add coverage of vertices to edges
# collapse graph
# compute edge coverage
# plot graph
a.fragmentate()
a.cover_edges()

print('Graph:')
print(a.graph_scheme)
print('Edges:')
print(a.edges)

a.edge_coverage()


a.plot('test', True)
a.collapse()
a.graph_scheme = a.collapsed_graph
a.edges = a.collapsed_edges
print(a.collapsed_graph)
print(a.edges)
a.plot('test_collapsed', True)


for v in a.graph_scheme.values():
    print(v[0])

for e in a.edges.values():
    print(e)
