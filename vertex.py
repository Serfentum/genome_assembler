from collections import defaultdict
from edge import Edge


class Vertex:
    # class for kmer vertices in dbgraph
    # Denotes kmer in a graph
    # Have
    # sequence - str
    # coverage - int
    # in degree -  dict with Vertex as a key and Edge as value
    # out degree
    # Maxima of in and out degrees is a number of letters in alphabet
    # function to update coverage of vertex
    def __init__(self, sequence, source=None, coverage=0):
        """
        Constructor
        :param sequence: str - sequence of kmer
        :param coverage: int - coverage of kmer
        :param source: Vertex - previous Vertex connected with this via Edge
        """
        self.sequence = sequence
        self.coverage = coverage

        # Perhaps redundant
        self.in_edges = defaultdict(list)
        if source:
            self.in_edges[source].append(Edge(sequence + source.sequence[-1]))
        self.out_edges = defaultdict(list)

    def __eq__(self, other):
        return other.__class__ == self.__class__ and self.sequence == other.sequence

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.sequence)

    def __str__(self):
        return str(self.__dict__)
