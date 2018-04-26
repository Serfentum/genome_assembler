from new.edge import Edge
from new.vertex import Vertex
from new.graph import Graph


# a = Graph('../test3', 2)
a = Graph('../wloop', 3)

a.fragmentate()
a.cover_edges()
a.edge_coverage()


a.plot('ttt', True, False)

a.collapse()
a.plot('ttr', False)
# a.graph_scheme = a.collapsed_graph
# a.edges = a.collapsed_edges
# print(a.collapsed_graph)
# print(a.edges)
a.plot('test_collapsed', include_seq=True, collapsed=True)

# print('Graph:')
# print(a.graph_scheme)
# print('Edges:')
# print(a.edges)
#
# for vs, (v, adj) in a.collapsed_graph.items():
#     print(v, adj, sep='\t')
# for e in a.edges.values():
#     print(e)



