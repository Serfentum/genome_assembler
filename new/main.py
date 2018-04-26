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
# for v in a.graph_scheme.values():
#     print(v[0])
# for e in a.edges.values():
#     print(e)



