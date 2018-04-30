# Genome Assembler
This is a genome assembler (GA) written in python. It can 
- build De Bruijn graph representation from fasta file with reads 
- compute coverages of vertices and edges
- collapse graph i.e. assemble genome
- filter low-quality vertices
- write assembly contigs to fasta
- draw graph at the following states: original, collapsed and stepwise - at each collapsing event
- create gif showing assembly

We have reworked collapse algorhitm and now it can process genomes ~10kb at short time. You can assemble genome by our assembler, but we would recommend to you use something more efficient e.g. [Spades](https://www.ncbi.nlm.nih.gov/pubmed/22506599).

The most unique feature of our assembler is ability to creates fiery gifs with assembly. This is commonly used with small genomes in educational purpose - to demonstrate how graph of assembly evolve.



## Prerequisites
1. python3.6 interpretator
2. biopython
3. graphviz
4. imageio


## Installation
To take a copy of our GA you should run this command in terminal

`git clone https://github.com/Serfentum/genome_assembler.git`

After this you can make use of GA with main.py as an entry point.

We are going to create cli and even gui interface


## Stucture
- src - source code
    - main - entry point
    - vertex - Vertex class
    - edge - Edge class
    - graph - Graph class
    - plots - plotting utilities
- examples - you can find awesome pictures and gifs with original and collapsed graph
    - hw - pictures and dots related to homework files
        - 1read - assembly of 1 read
        - 1virus - assembly of 1 virus
- additional_tests - some other tests

## Usage
For now you should go into main file and look at `launch` function structure. It has calls to functions in appropriate order to create a graph from reads, collapse it, remove low-covered vertices, plot and create gif.

There are several parameters to launch functions:
- filepath - path to input fasta file
- k - size of kmer
- include_seq - whether to include sequences of vertex and edges in plot
- threshold - level of noise-tolerance - the more threshold the more coverage value will be needed for vertex to stay at graph; threshold = 1 will delete all vertices with coverage lower than mean
- format - format of plots, pdf, svg and png are appropriate; for gif only png is acceptable one
- output - path to output directory
- show - whether to show you all plots in a process of genome assembling
- pause - time between each plot in show and gif
- fix_steps - whether to create plot at each assembly step
- make_gif - whether to create gif

## Versioning
0.3.0


## Authors
* Alexander Ilin

## Acknowledgements
I'd like to thank Nastya, Jenya, Anton and Sasha.

## License
BioInformatics license
