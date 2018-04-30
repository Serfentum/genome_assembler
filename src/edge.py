import numpy as np


class Edge:
    k = None

    def __init__(self, source, dest):
        self.source = source
        self.dest = dest
        self.sequence = source.sequence + dest.sequence[-1]
        self.coverage = 0
        self.coverages = 0
        self.parts = 2

    def ini_cover(self):
        """
        Place source and dest coverages into coverages array
        :return:
        """
        self.coverages = self.source.coverage + self.dest.coverage
        self.coverage = self.coverages / self.parts

    def compute_coverage(self):
        """
        Compute coverage of edge as a mean of including vertex coverages
        :return:
        """
        self.coverage = self.coverages / self.parts

    def extend(self, second_edge):
        """
        Modify edge - add last coverage from 2nd edge to coverages, extend sequence and update dest
        :param second_edge: Edge - 2nd edge to the given to form extended edge
        :return:
        """
        self.coverages += second_edge.coverages - second_edge.source.coverage
        self.sequence += second_edge.sequence[Edge.k:]
        # self.sequence += second_edge.sequence[-(len(second_edge.coverages) - 1):]
        self.dest = second_edge.dest
        self.parts += second_edge.parts - 1
        # print(f'{self.sequence}\t{self.coverages}\t{self.parts}')

    def __eq__(self, other):
        return other.__class__ == self.__class__ and self.sequence == other.sequence

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.sequence)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)
