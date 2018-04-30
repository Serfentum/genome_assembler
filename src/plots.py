import time
from graphviz import Digraph


def plot(graph, filename, include_seq, format, show):
    """
    Plot assembly graph
    :param graph: Graph - object of graph
    :param filename: str - path to output
    :param include_seq: boolean - whether to include sequences in vertices and edges
    :param format: str - format of slides
    :return:
    """
    # Graph creation and parameters setting
    dot = Digraph(comment='Assembly', format=format)
    dot.graph_attr.update(size='7.75,10.25!', ratio='fill')

    # Choose function with appropriate labeling according to full or short plot
    if include_seq:
        dot = plot_full(graph, dot)
    else:
        dot = plot_wo_seq(graph, dot)
    # Save pdf and perhaps show it
    if show:
        dot.render(filename, view=True)
    else:
        dot.render(filename)


def slideshow(graph, pause, format, show, output):
    """
    Plot current state of graph, nice to use in loops to visualize
    :param graph: Graph - object of graph
    :param pause: float - seconds to display picture
    :param format: str - format of slides
    :return:
    """
    # Graph creation and parameters setting
    dot = Digraph(comment='Assembly', format=format)
    dot.graph_attr.update(size='7.75,10.25!', ratio='fill')

    # Plot graph, display it and wait pause time if should
    plot_full(graph, dot)
    if show:
        dot.render(f'{output}/picts/{time.time()}', view=True)
        time.sleep(pause)
    else:
        dot.render(f'{output}/picts/{time.time()}')


def plot_wo_seq(graph, dot):
    """
    Create content of dot file with short graph description without sequence labels
    :param graph: Graph - object of graph
    :param dot: Digraph - object from graphviz to describe graph
    :return: Digraph - fulfilled with nodes and edges
    """
    # Iterate over vertices and edges, add them to graphviz graph
    # There will be coverages of vertices, edges and their length
    for vs, (v, out) in graph.graph_scheme.items():
        dot.node(vs, label=f'{v.coverage}')
        for dv in out:
            for e in graph.edges[(vs, dv)]:
                dot.edge(vs, dv, label=f'{graph.k + len(e.coverages) - 1} {e.coverage}')
    return dot


def plot_full(graph, dot):
    """
    Create content of dot file with full graph description including sequence labels
    :param graph: Graph - object of graph
    :param dot: Digraph - object from graphviz to describe graph
    :return: Digraph - fulfilled with nodes and edges
    """
    # Iterate over vertices and edges, add them to graphviz graph
    # There will be coverages of vertices, edges and their length and corresponding sequences
    for vs, (v, out) in graph.graph_scheme.items():
        dot.node(vs, label=f'{vs} {v.coverage}')
        for dv in set(out):
            for e in set(graph.edges[(vs, dv)]):
                dot.edge(vs, dv, label=f'{e.sequence} {graph.k + e.parts - 1} {e.coverage}')
    return dot