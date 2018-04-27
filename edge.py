import numpy as np


class Edge:
    def __init__(self, source, dest):
        self.source = source
        self.dest = dest
        self.sequence = source.sequence + dest.sequence[-1]
        self.coverage = 0
        self.coverages = np.array(())

    def ini_cover(self):
        """
        Place source and dest coverages into coverages array
        :return:
        """
        self.coverages = np.append(self.coverages, (self.source.coverage, self.dest.coverage))

    def compute_coverage(self):
        """
        Compute coverage of edge as a mean of including vertex coverages
        :return:
        """
        self.coverage = self.coverages.mean()

    def extend(self, second_edge):
        """
        Modify edge - add last coverage from 2nd edge to coverages, extend sequence and update dest
        :param second_edge: Edge - 2nd edge to the given to form extended edge
        :return:
        """
        self.coverages = np.append(self.coverages, (second_edge.coverages[1:]))
        self.sequence += second_edge.sequence[-(len(second_edge.coverages) - 1):]
        self.dest = second_edge.dest

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
