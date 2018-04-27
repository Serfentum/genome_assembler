from graph import Graph


def launch(filepath, k, include_seq=True, threshold=0.3, format='pdf', output='out', show=False, pause=1.2):
    """
    Launch function
    :return:
    """
    # Load file
    a = Graph(filepath, k)

    # Populate graph
    a.fragmentate()
    # Initial edges covering
    a.cover_edges()
    a.edge_coverage()

    # Make plot of full graph
    a.plot(f'{output}_original', include_seq, format)

    # Collapse graph and removing low covered vertices
    a.collapse(show, pause, format)
    a.remove_outliers(threshold)
    # Compute edge coverage
    a.edge_coverage()

    # Create plot of collapsed graph
    a.plot(f'{output}collapsed', include_seq, format)



# Sets
# a = Graph('../linear', 21)
# a = Graph('/home/arleg/Downloads/hw3_dataset.fasta', 55)
# a = Graph('/home/arleg/Downloads/hw3_dataset(1).fasta', 3)
# a = launch('../wloop', 3, show=True, format='svg', pause=0.5)
a = launch('../test3', 3)
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



