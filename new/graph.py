from graphviz import Digraph
import logging
import time
from Bio import SeqIO
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
        self.graph_scheme = {}
        self.edges = {}
        self.collapsed_graph = {}
        self.collapsed_edges = {}

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
                dot.edge(vs, dv,
                         label=f'{self.k + len(edges[(vs, dv)].coverages) - 1} '
                               f'{edges[(vs, dv)].coverage}')
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
            for dv in out:
                dot.edge(vs, dv,
                         label=f'{edges[(vs, dv)].sequence} '
                               f'{self.k + len(edges[(vs, dv)].coverages) - 1} '
                               f'{edges[(vs, dv)].coverage}')
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
        logger.debug('Collapsing\n')
        while True:
            self.collapse_iteration(show, pause, format)
            logger.debug('--> Number of vertices in graph is %s', len(self.graph_scheme))
            if self.graph_scheme == self.collapsed_graph and self.edges == self.collapsed_edges:
                break
            self.graph_scheme = self.collapsed_graph
            self.edges = self.collapsed_edges
            self.collapsed_graph = {}
            self.collapsed_edges = {}
            if len(self.graph_scheme) == 1:
                break

    def collapse_iteration(self, show, pause, format):
        """
        Make 1 collapsing iteration
        Populate self.collapsed_graph and self.collapsed_edges
        :param show: boolean - whether to visualize each step
        :param pause: float - time of pauses between slides
        :param format: str - format of slides
        :return:
        """
        obsolete = set()
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

            for inter in adj_vertex:
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
                        continue
                    logger.debug('\t\t\t\tExtend edge with source, inter, dest %s %s %s', source.sequence, inter.sequence, dest.sequence)
                    self.extend_one_edge(source, inter, dest)
                logger.debug('\t\tAdd %s to obsolete', inter.sequence)
                obsolete.add(inter.sequence)
                if show:
                    self.slideshow(pause, True, format)

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
        """
        Write information about node coverage into edges
        :return:
        """
        for edge in self.edges.values():
            edge.ini_cover()

    def edge_coverage(self):
        """
        Compute coverage of edges
        :return:
        """
        for edge in self.edges.values():
            edge.compute_coverage()

    def remove_outliers(self, threshold):
        """
        Remove from graph vertices and edges with low coverage
        :param threshold: float - all vertices with coverage < threshold * mean(coverage) will be removed
        :return:
        """
        cutoff = threshold * self.mean_coverage()
        # Select vertices with coverage less than cutoff
        # Condition about out degree == 0 should be added if you want to delete just leaves
        obsolete_vertices = (i[0] for i in (filter(lambda v: getattr(v[0], 'coverage') < cutoff, self.graph_scheme.values())))

        # Delete vertices with low coverage and their edges
        for v in list(obsolete_vertices):
            del self.graph_scheme[v.sequence]
            # Select all edges where obsolete vertex present
            obsolete_edges = [e for e in self.edges if e[0] == v.sequence or e[1] == v.sequence]
            # Delete obsolete edge from edges and delete links in graph
            for e in obsolete_edges:
                del self.edges[e]
                try:
                    self.graph_scheme[e[0]][1].remove(v.sequence)
                    self.graph_scheme[e[1]][1].remove(v.sequence)
                except:
                    pass


    def mean_coverage(self):
        """
        Compute mean vertex coverage
        :return: float - mean vertex coverage
        """
        return sum(getattr(v[0], 'coverage') for v in self.graph_scheme.values()) / len(self.graph_scheme)

