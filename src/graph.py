from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from graphviz import Digraph
import logging
import time
import imageio
from edge import Edge
from vertex import Vertex


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

    def plot(self, filename, include_seq, format, collapsed=False):
        """
        Plot assembly graph
        :param filename: str - path to output
        :param include_seq: boolean - whether to include sequences in vertices and edges
        :param format: str - format of slides
        :param collapsed: boolean - whether collapsed graph to plot #TODO remove this parameter
        :return:
        """
        dot = Digraph(comment='Assembly', format=format)

        # Choose function with appropriate labeling according to full or short plot
        if include_seq:
            dot = self.plot_full(dot, collapsed)
        else:
            dot = self.plot_wo_seq(dot, collapsed)
        # Save pdf
        dot.render(filename, view=True)

    def slideshow(self, pause, collapsed, format):
        """
        Plot current state of graph, nice to use in loops to visualize
        :param pause: float - seconds to display picture
        :param format: str - format of slides
        :param collapsed: boolean #TODO remove obsolete parameter
        :return:
        """
        dot = Digraph(comment='Assembly', format=format)

        # Plot graph, display it and wait pause time
        self.plot_full(dot, collapsed)
        dot.render(view=True)
        time.sleep(pause)

    def plot_wo_seq(self, dot, collapsed=False):
        """
        Create content of dot file with short graph description without sequence labels
        :param dot: Digraph - object from graphviz to describe graph
        :param collapsed: boolean - whether collapsed graph to plot
        :return: Digraph - fulfilled with nodes and edges
        """
        # Choose appropriate graph for plotting - full or collapsed
        graph, edges = self.solve_graph_type(collapsed)

        # Iterate over vertices and edges, add them to graphviz graph
        # There will be coverages of vertices, edges and their length
        for vs, (v, out) in graph:
            dot.node(vs, label=f'{v.coverage}')
            for dv in out:
                for e in edges[(vs, dv)]:
                    dot.edge(vs, dv, label=f'{self.k + len(e.coverages) - 1} {e.coverage}')
        return dot

    def plot_full(self, dot, collapsed=False):
        """
        Create content of dot file with full graph description including sequence labels
        :param dot: Digraph - object from graphviz to describe graph
        :param collapsed: boolean - whether collapsed graph to plot
        :return: Digraph - fulfilled with nodes and edges
        """
        # Choose appropriate graph for plotting - full or collapsed
        graph, edges = self.solve_graph_type(collapsed)

        # Iterate over vertices and edges, add them to graphviz graph
        # There will be coverages of vertices, edges and their length and corresponding sequences
        for vs, (v, out) in graph:
            dot.node(vs, label=f'{vs} {v.coverage}')
            for dv in set(out):
                for e in set(edges[(vs, dv)]):
                    dot.edge(vs, dv, label=f'{e.sequence} {self.k + e.parts - 1} {e.coverage}')
        return dot

    def solve_graph_type(self, collapsed):
        """
        Choose appropriate graph for plotting - full or collapsed
        :param collapsed: boolean - whether collapsed graph to plot
        :return: dict, dict - tuple of appropriate graph and edges dicts
        """
        if collapsed:
            graph = self.collapsed_graph.items()
            edges = self.collapsed_edges
        else:
            graph = self.graph_scheme.items()
            edges = self.edges
        return graph, edges

    def collapse(self, show, pause, format):
        """
        Collapse graph, after this self.collapsed_graph will be in self.graph_scheme
        :param show: boolean - whether to visualize each step
        :param pause: float - time of pauses between slides
        :param format: str - format of slides
        :return:
        """
        # cycle indeterminate
        # collapse iteration
        # reassignments:
        # Find all collapsable vertices i.e. those which has in and out degrees equal to 1
        collapsable = list(filter(lambda vertex: vertex.in_degree == 1 and vertex.out_degree == 1, (x[0] for x in self.graph_scheme.values())))
        # Display uncollapsed graph
        if show:
            self.slideshow(pause, False, format)

        # Iterate over collapsed vertices
        # reassign edges from previous vertices to next in self.edges and self.graph_scheme
        # update in and out vertices in vertices from both sides
        # delete collapsed from self.edges and self.graph_scheme
        logging.debug('Start collapsing:')
        for collapsing in collapsable:
            # print(collapsing, sep='\t')
            # Relink vertices in graph_scheme
            self.graph_scheme[collapsing.ins[0]][1].remove(collapsing.sequence)
            self.graph_scheme[collapsing.ins[0]][1].append(collapsing.outs[0])

            # Vertex in and out neighbours update
            self.graph_scheme[collapsing.ins[0]][0].outs.remove(collapsing.sequence)
            self.graph_scheme[collapsing.ins[0]][0].outs.append(collapsing.outs[0])
            self.graph_scheme[collapsing.outs[0]][0].ins.remove(collapsing.sequence)
            self.graph_scheme[collapsing.outs[0]][0].ins.append(collapsing.ins[0])

            # Edge update
            first_edge = self.edges[(collapsing.ins[0], collapsing.sequence)][0]
            second_edge = self.edges[(collapsing.sequence, collapsing.outs[0])][0]
            first_edge.extend(second_edge)
            self.edges[(collapsing.ins[0], collapsing.outs[0])].append(first_edge)

            # Collapsed vertex and edge removing
            del self.graph_scheme[collapsing.sequence]
            del self.edges[(collapsing.sequence, collapsing.outs[0])]
            del self.edges[(collapsing.ins[0], collapsing.sequence)]

            # Show image of graph
            if show:
                self.slideshow(pause, False, format)

    def extend_one_edge(self, source, inter, dest):
        """
        Extend (source, inter) to (source, dest)
        :param source:
        :param inter:
        :param dest:
        :return:
        """
        # Add source vertex to collapsed graph
        if source.sequence not in self.collapsed_graph:
            self.collapsed_graph[source.sequence] = source, []

        # Add link between 1st and 3rd vertices
        self.collapsed_graph[source.sequence][1].append(dest.sequence)
        # Add edge to collapsed edges
        self.modify_edge(source, inter, dest)

        # If inter vertex was in collapsed graph delete it
        if inter.sequence in self.collapsed_graph:
            del self.collapsed_graph[inter.sequence]
            logger.debug('%s', self.collapsed_graph.keys())

    def modify_edge(self, source, inter, dest):
        """
        Assign edge between source and dest vertices in collapsed graph, update edge's sequence, coverages and dest
        :param source: Vertex - 1st layer vertex
        :param inter: Vertex - 2nd layer vertex, next to source
        :param dest: Vertex - 3rd layer vertex, next to inter
        :return:
        """
        e = self.edges[(source.sequence, inter.sequence)][0]
        self.collapsed_edges[(source.sequence, dest.sequence)].append(e)
        e.extend(self.edges[(inter.sequence, dest.sequence)][0])

    def add_unchanged(self, source, inter=None):
        """
        Add source vertex or source vertex with inter and their edge into new graph intact
        :param source: Vertex - first vertex
        :param inter: Vertex - second vertex
        :return:
        """
        if source.sequence not in self.collapsed_graph:
            self.collapsed_graph[source.sequence] = source, []
        if inter:
            if inter.sequence not in self.collapsed_graph:
                self.collapsed_graph[inter.sequence] = inter, []
            self.collapsed_graph[source.sequence][1].append(inter.sequence)
            self.collapsed_edges[(source.sequence, inter.sequence)] = self.edges[(source.sequence, inter.sequence)]



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
            self.edges[(source, destination)].append(Edge(self.graph_scheme[source][0], self.graph_scheme[destination][0]))

    def cover_edges(self):
        """
        Write information about node coverage into edges
        :return:
        """
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
        cutoff = threshold * self.mean_coverage()
        # Select vertices with coverage less than cutoff
        # Condition about out degree == 0 should be added if you want to delete just leaves
        # Deleting below will be easier in this case
        obsolete_vertices = list(inter[0] for inter
                             in (filter(lambda vertex: getattr(vertex[0], 'coverage') < cutoff,
                                        self.graph_scheme.values())))
        # obsolete_edges = set()
        # Delete vertices with low coverage and their edges
        # Delete links in graph
        # Update adjacent vertex in and out degrees
        for vertex in list(obsolete_vertices):
            for previous in vertex.ins:
                self.graph_scheme[previous][0].outs.remove(vertex.sequence)
                self.graph_scheme[previous][1].remove(vertex.sequence)
                self.graph_scheme[previous][0].out_degree -= 1
                # obsolete_edges.add()
                del self.edges[(previous, vertex.sequence)]
            # print(obsolete_edges)
            for following in vertex.outs:
                # try:
                self.graph_scheme[following][0].ins.remove(vertex.sequence)
                self.graph_scheme[following][0].in_degree -= 1
                # obsolete_edges.add((vertex.sequence, following))
                del self.edges[(vertex.sequence, following)]
                # except:
                #     pass
                # self.graph_scheme[following][1].remove(vertex.sequence)

            del self.graph_scheme[vertex.sequence]

            # Select all edges where obsolete vertex present
            # obsolete_edges = [e for e in self.edges if e[0] == v.sequence or e[1] == v.sequence]
            # Delete obsolete edge from edges and delete links in graph
            # for edge in obsolete_edges:
                # try:

                    # self.graph_scheme[edge[0]][1].remove(vertex.sequence)
                    # self.graph_scheme[edge[1]][1].remove(vertex.sequence)
                # except:
                #     pass

    def mean_coverage(self):
        """
        Compute mean vertex coverage
        :return: float - mean vertex coverage
        """
        return sum(getattr(v[0], 'coverage') for v in self.graph_scheme.values()) / len(self.graph_scheme)

    def extract(self, output):
        # Simple variant - dump all edges
        edges = [SeqRecord(Seq(getattr(e, 'sequence')), description='', id='edge_from_ga')
                 for edge in self.edges.values() for e in edge]
        SeqIO.write(edges, f'{output}.fasta', 'fasta')
