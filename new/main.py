from new.edge import Edge
from new.vertex import Vertex
from new.graph import Graph


# a = Graph('../linear', 21)
# a = Graph('/home/arleg/Downloads/hw3_dataset.fasta', 55)
# a = Graph('/home/arleg/Downloads/hw3_dataset(1).fasta', 3)
# a = Graph('../wloop', 3)
a = Graph('../test4', 3)

a.fragmentate()
a.cover_edges()
a.edge_coverage()


a.plot('linear_test2', True, False, format='png')

a.collapse()
# a.plot('linear_test_col', True, False)
# a.plot('ttr', False, False)
# a.graph_scheme = a.collapsed_graph
# a.edges = a.collapsed_edges
# print(a.collapsed_graph)
# print(a.edges)
# a.plot('test_collapsed', include_seq=True, collapsed=False)

# print('Graph:')
# print(a.graph_scheme)
# print('Edges:')
# print(a.edges)
#
# for vs, (v, adj) in a.collapsed_graph.items():
#     print(v, adj, sep='\t')
# for e in a.edges.values():
#     print(e)



