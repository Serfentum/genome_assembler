import numpy as np


class Edge:
    def __init__(self, source, dest):
        self.source = source
        self.dest = dest
        self.sequence = source.sequence + dest.sequence[-1]
        self.coverage = 0
        self.coverages = np.array(())

    def ini_cover(self):
        self.coverages = np.append(self.coverages, (self.source.coverage, self.dest.coverage))

    def compute_coverage(self):
        self.coverage = self.coverages.mean()

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
