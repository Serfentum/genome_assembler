from graph import Graph


def launch(filepath, k, include_seq=True, threshold=0.3, format='pdf', output='out', show=False, pause=1.2):
    """
    Launch function to create graph representation of assembly, collapse it, remove low covered vertices
    :param filepath: str - path to fasta file
    :param k: int - kmer size, recommended to be odd
    :param include_seq: boolean - whether to display vertex and edge sequences
    :param threshold: float - coefficient which determines threshold for coverage for vertex removing
    :param format: str - format of output images - pdf, svg, png
    :param output: str - file name of output, 'original' and 'collapsed' will be appended to it
    :param show: boolean - whether to display fiery slide show
    :param pause: float - pause between each slide
    :return:
    """
    # Load file
    a = Graph(filepath, k)

    # Populate graph and cover edges
    a.fragmentate()
    a.cover_edges()
    a.edge_coverage()

    # Make plot of full graph
    # a.plot(f'{output}/original', include_seq, format)

    # Collapse graph and removing low covered vertices
    print(len(a.graph_scheme))
    for i in range(3):
        a.collapse(show, pause, format)
        a.remove_outliers(threshold)
        print(len(a.graph_scheme))
        print(len(list(filter(lambda vertex: vertex.in_degree == 1 and vertex.out_degree == 1, (x[0] for x in a.graph_scheme.values())))))

    # Compute edge coverage
    a.edge_coverage()
    a.extract(f'{output}/fasta')

    # Create plot of collapsed graph
    a.plot(f'{output}/collapsed', include_seq, format)



# Sets
# launch('../examples/linear', 21, show=False)
# launch('../examples/linear', 21, threshold=0.5, show=True, format='pdf', pause=0.4, output='/home/arleg/genome_assembler/out')
launch('/home/arleg/Downloads/hw3_dataset(1).fasta', 55, show=False, output='/home/arleg/genome_assembler/test/out/05', threshold=0.5)
# a = launch('examples/snp', 3, show=True, threshold=0.1, format='pdf', pause=1)
# launch('../additional_tests/test3', 3, threshold=0.7, output='/home/arleg/genome_assembler/test/out', show=True)
# a = Graph('/home/arleg/Downloads/s_6_first100000.fastq', 55)


# a.plot('Nastya_test', True, False, format='pdf')
# a.collapse()
# a.edge_coverage()
# a.remove_outliers(0.7)
# a.plot('Nastya_test_col', True, False, format='pdf')
# a.plot('linear_test_col', True, False)


# for vs, (v, adj) in a.collapsed_graph.items():
#     print(v, adj, sep='\t')
# for e in a.edges.values():
#     print(e)
#     print(e.coverages.mean(), e.coverage)



