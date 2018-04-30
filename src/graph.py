from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
from edge import Edge
from vertex import Vertex
from plots import slideshow


# For logging purpose. It will create file graph.log
# To turn logging on    - change logging level to DEBUG
# To turn logging off   - change logging level to something higher e.g. INFO
logging.basicConfig(filename='graph.log',
                    filemode='w',
                    format='%(message)s',
                    level=logging.INFO)
logger = logging.getLogger(__name__)


class Graph:
    def __init__(self, path, k):
        self.reads = SeqIO.parse(path, 'fasta')
        self.k = k
        Edge.k = k
        self.graph_scheme = {}
        self.vertices = set()
        self.edges = defaultdict(list)
        self.collapsed_graph = {}
        self.collapsed_edges = defaultdict(list)

    def collapse_filter(self, threshold, fix_steps, show, pause, format, output):
        """
        Collapse graph and remove low covered vertices after this collapse it again - to most possible collapsed state
        :param threshold: float - all vertices with coverage < threshold * mean(coverage) will be removed
        :param show: boolean - whether to visualize each step
        :param pause: float - time of pauses between slides
        :param format: str - format of slides
        :return:
        """
        # Initial collapsing
        logging.debug('Original number of vertices is %s', len(self.graph_scheme))
        self.collapse(fix_steps, show, pause, format, output)
        logging.debug('Number of vertices after 1st collapsing - %s', len(self.graph_scheme))

        # Removing low covered
        self.remove_outliers(threshold)
        logging.debug('Number of vertices after removing low-covered - %s', len(self.graph_scheme))

        # Collapsing
        self.collapse(fix_steps, show, pause, format, output)
        logging.debug('Number of vertices after next collapsing - %s', len(self.graph_scheme))

    def collapse(self, fix_steps, show, pause, format, output):
        """
        Collapse graph, after this self.collapsed_graph will be in self.graph_scheme
        :param fix_steps: boolean - whether to visualize each step
        :param pause: float - time of pauses between slides
        :param format: str - format of slides
        :return:
        """
        # Find all collapsable vertices i.e. those which has in and out degrees equal to 1
        collapsable = list(filter(lambda vertex: vertex.in_degree == 1 and vertex.out_degree == 1, (x[0] for x in self.graph_scheme.values())))
        logging.debug('Collapsable vertex are:\n%s', collapsable)
        # Display uncollapsed graph
        if fix_steps:
            slideshow(self, pause=pause, format=format, show=show, output=output)

        # Iterate over collapsed vertices
        # reassign edges from previous vertices to next in self.edges and self.graph_scheme
        # update in and out vertices in vertices from both sides
        # delete collapsed from self.edges and self.graph_scheme
        # Show state of graph
        logging.debug('Start collapsing:')
        if fix_steps:
            for collapsing in collapsable:
                self.collapse_vertex(collapsing)
                slideshow(self, pause=pause, format=format, show=show, output=output)
        else:
            for collapsing in collapsable:
                self.collapse_vertex(collapsing)

    def collapse_vertex(self, collapsing):
        """
        Collapse 1 vertex, extend it source to destination, update structures
        :param collapsing: Vertex - collapsed vertex
        :return:
        """
        # Update vertices, edges and delete collapsed vertex and edges
        self.update_vertex(collapsing)
        self.update_edge(collapsing)
        self.clean(collapsing)

    def update_vertex(self, collapsing):
        """
        Update vertex links in graph and Vertex ins and outs
        :param collapsing: Vertex - collapsed vertex
        :return:
        """
        # Relink vertices in graph_scheme
        self.graph_scheme[collapsing.ins[0]][1].remove(collapsing.sequence)
        self.graph_scheme[collapsing.ins[0]][1].append(collapsing.outs[0])

        # Vertex in and out neighbours update
        self.graph_scheme[collapsing.ins[0]][0].outs.remove(collapsing.sequence)
        self.graph_scheme[collapsing.ins[0]][0].outs.append(collapsing.outs[0])
        self.graph_scheme[collapsing.outs[0]][0].ins.remove(collapsing.sequence)
        self.graph_scheme[collapsing.outs[0]][0].ins.append(collapsing.ins[0])

    def update_edge(self, collapsing):
        """
        Unite 2 edges together, update coverage and self.edges
        :param collapsing: Vertex - collapsed vertex
        :return:
        """
        # Take necessary edges
        first_edge = self.edges[(collapsing.ins[0], collapsing.sequence)][0]
        second_edge = self.edges[(collapsing.sequence, collapsing.outs[0])][0]

        # Modify edge and place it into edge dictionary
        first_edge.extend(second_edge)
        self.edges[(collapsing.ins[0], collapsing.outs[0])].append(first_edge)

    def clean(self, collapsing):
        """
        Delete traces of collapsed vertex from graph and edges dictionaries
        :param collapsing: Vertex - collapsed vertex
        :return:
        """
        # Collapsed vertex and edge removing
        del self.graph_scheme[collapsing.sequence]
        del self.edges[(collapsing.sequence, collapsing.outs[0])]
        del self.edges[(collapsing.ins[0], collapsing.sequence)]

    def fragmentate(self):
        """
        Cleave reads into kmers and add them to graph
        :return:
        """
        # Cleave all reads to kmers
        for read in self.reads:
            self.add_read(str(read.seq).upper())
            self.add_read(str(read.reverse_complement().seq).upper())

    def add_read(self, read):
        """
        Add kmers from read to graph
        :param read: str - read sequence
        :return:
        """
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

                    self.graph_scheme[source][0].outs.append(destination)
                    self.graph_scheme[destination][0].ins.append(source)

                    edge = Edge(self.graph_scheme[source][0], self.graph_scheme[destination][0])
                    self.edges[(source, destination)].append(edge)
                continue

            # Create absent vertex
            elif source in self.graph_scheme and destination not in self.graph_scheme:
                self.graph_scheme[destination] = (Vertex(destination), [])
                self.graph_scheme[source][0].outs.append(destination)
                self.graph_scheme[destination][0].ins.append(source)
            elif source not in self.graph_scheme and destination in self.graph_scheme:
                self.graph_scheme[source] = (Vertex(source), [])
                self.graph_scheme[destination][0].ins.append(source)
                self.graph_scheme[source][0].outs.append(destination)
            else:
                self.graph_scheme[source] = (Vertex(source), [])
                self.graph_scheme[destination] = (Vertex(destination), [])
                self.graph_scheme[source][0].outs.append(destination)
                self.graph_scheme[destination][0].ins.append(source)

            # ss =
            # self.vertices.add(Vertex(source))
            # self.vertices.add(Vertex(destination))

            # Increase degrees of vertices
            self.graph_scheme[source][0].increase_out_degree()
            self.graph_scheme[destination][0].increase_in_degree()
            # Increase vertex coverage
            self.graph_scheme[source][0].increase_coverage()
            self.graph_scheme[destination][0].increase_coverage()
            # Add link in graph edge to edges set
            self.graph_scheme[source][1].append(destination)
            self.edges[(source, destination)].append(Edge(self.graph_scheme[source][0],
                                                          self.graph_scheme[destination][0]))

    def cover_edges(self):
        """
        Write information about node coverage into edges
        :return:
        """
        # Assign coverage to each edge as a ratio of including vertex coverages to their number
        for edge in self.edges.values():
            for e in edge:
                e.ini_cover()

    def edge_coverage(self):
        """
        Compute coverage of edges
        :return:
        """
        for edge in self.edges.values():
            for e in edge:
                e.compute_coverage()

    def remove_outliers(self, threshold):
        """
        Remove from graph vertices and edges with low coverage
        :param threshold: float - all vertices with coverage < threshold * mean(coverage) will be removed
        :return:
        """
        # Threshold for determining low covered vertices
        cutoff = threshold * self.mean_coverage()

        # Select vertices with coverage less than cutoff
        # Conditioning on out_degree == 0 should be added if you want to delete just leaves
        # Deleting below will be easier in this case
        obsolete_vertices = list(inter[0] for inter in (filter(lambda vertex: getattr(vertex[0], 'coverage') < cutoff,
                                                               self.graph_scheme.values())))
        logging.debug('Low covered vertex are:\n%s', obsolete_vertices)
        # Iterate over vertices with low coverage, for vertices to which these are adjacent and for adjacent vertices
        # Delete links in Vertex, links in graph, update in and out degrees, delete their edges and vertices itself
        for vertex in obsolete_vertices:
            for previous in vertex.ins:
                self.graph_scheme[previous][0].outs.remove(vertex.sequence)
                self.graph_scheme[previous][1].remove(vertex.sequence)
                self.graph_scheme[previous][0].out_degree -= 1
                del self.edges[(previous, vertex.sequence)]

            for following in vertex.outs:
                self.graph_scheme[following][0].ins.remove(vertex.sequence)
                self.graph_scheme[following][0].in_degree -= 1
                del self.edges[(vertex.sequence, following)]

            del self.graph_scheme[vertex.sequence]

    def mean_coverage(self):
        """
        Compute mean vertex coverage
        :return: float - mean vertex coverage
        """
        return sum(getattr(v[0], 'coverage') for v in self.graph_scheme.values()) / len(self.graph_scheme)

    def extract(self, output):
        """
        Dump all edges in a fasta
        :param output: str - name of output file
        :return:
        """
        # Take all edges
        edges = [SeqRecord(Seq(getattr(e, 'sequence')), description='', id='edge_from_ga')
                 for edge in self.edges.values() for e in edge]
        # Write edges to fasta file
        SeqIO.write(edges, f'{output}.fasta', 'fasta')
