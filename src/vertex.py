class Vertex:
    def __init__(self, sequence):
        self.sequence = sequence
        self.coverage = 0
        # seqs
        self.ins = []
        self.outs = []
        self.in_degree = 0
        self.out_degree = 0

    def increase_coverage(self, n=1):
        self.coverage += n

    def increase_in_degree(self):
        self.in_degree += 1

    def increase_out_degree(self):
        self.out_degree += 1

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
