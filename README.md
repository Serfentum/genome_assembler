# Genome Assembler
This is a genome assembler (GA) written in python. It is overengineered so much, that can't handle large enough fasts now. Probably due to overly redundant structures and non-efficient collapsing algorhitm. 

But you still can construct graphs, compute coverages, collapse graph, prune low-covered elements and look at fiery slide shows for small files!

Feel free to advise how to lighten graph structure.


## Prerequisites
1. python3.6 interpretator
2. biopython
3. graphviz


## Installation
To take a copy of our GA you should run this command in terminal

`git clone https://github.com/Serfentum/genome_assembler.git`

After this you can make use of GA with main.py as an entry point.


## Stucture
- src - source code
    - main - entry point
    - vertex - Vertex class
    - edge - Edge class
    - graph - Graph class
- examples - you can find awesome pictures and gifs with original and collapsed graph
    - hw - pictures and dots related to homework files
- additional_tests - some other tests


## Versioning
0.2.0


## Authors
* Alexander Ilin

## License
BioInformatics license
