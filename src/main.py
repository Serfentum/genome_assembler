from graph import Graph
from plots import plot
from animator import animate


def launch(filepath,
           k=55,
           include_seq=True,
           threshold=0.3,
           format='pdf',
           output='out',
           show=False,
           pause=1.2,
           fix_steps=False,
           make_gif=False):
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
    :param fix_steps: boolean - whether to create image of graph state at every stage
    :return:
    """
    # Load file
    a = Graph(filepath, k)

    # Populate graph and cover edges
    a.fragmentate()
    a.cover_edges()
    a.edge_coverage()

    # Collapse graph and remove low covered vertices
    a.collapse_filter(threshold, fix_steps=fix_steps, show=show, pause=pause, format=format, output=output)

    # Compute edge coverage
    a.edge_coverage()
    a.extract(f'{output}/out')

    # Create plot of collapsed graph if it wasn't created
    if not fix_steps:
        plot(a, f'{output}/picts/collapsed', include_seq, format, show)

    # Create animation from obtained images
    if make_gif:
        animate(f'{output}/picts', format, pause)



# Sets
# launch('../examples/variation', 3, fix_steps=True, show=True, threshold=0.5, output='var', format='png', make_gif=True)
# launch('../examples/linear', 21, threshold=0.5, show=True, format='pdf', pause=0.4, output='/home/arleg/genome_assembler/out')
# launch('/home/arleg/Downloads/hw3_dataset(1).fasta', 55, show=False, output='/home/arleg/genome_assembler/test/out/05', threshold=0.5)
# a = launch('examples/snp', 3, show=True, threshold=0.1, format='pdf', pause=1)
# launch('../additional_tests/test3', 3, threshold=0.7, output='/home/arleg/genome_assembler/test/out', show=True)
# a = Graph('/home/arleg/Downloads/s_6_first100000.fastq', 55)


launch('/home/arleg/genome_assembler/examples/hw/hw3_dataset.fasta',
           k=55,
           include_seq=True,
           threshold=0.3,
           format='png',
           output='/home/arleg/genome_assembler/examples/hw/1read',
           show=False,
           pause=0.3,
           fix_steps=True,
           make_gif=True)


