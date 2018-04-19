class Edge:
    # Have sequence (seq of 1st kmer and last nucleotide of 2nd),
    # coverage (1 by default, should be updated as a mean of 2 vertices)
    # number of kmers which comprise edge (default is 2)
    # function to update coverage of edge

    def __init__(self, sequence):
        self.sequence = sequence
        self.coverage = 0
        self.size = 2

    def increase_coverage(self, coverage=1):
        self.coverage += coverage

    def __eq__(self, other):
        return other.__class__ == self.__class__ and self.sequence == other.sequence

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.sequence)

    def __str__(self):
        return str(self.__dict__)
