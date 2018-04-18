class Edge:
    # Have sequence (seq of 1st kmer and last nucleotide of 2nd),
    # coverage (1 by default, should be updated as a mean of 2 vertices)
    # number of kmers which comprise edge (default is 2)
    # function to update coverage of edge

    def __init__(self, sequence):
        self.sequence = sequence
        self.coverage = 0
        self.size = 2

