<h1>The Problem</h1>

Experimental evidence has demonstrated that at least one large portion of the C. reinhardtii genome is misaligned in the current reference genome, with others suspected. The problem has persisted across multiple drafts of the genome.

<h1>The Solution (?)</h1>
A large library of BAC and fosmid paired-end reads aligned to the reference genome is available on phytozome. This dataset was downloaded and processed to remove the BAC data. Since a very high percentage of fosmid paired ends are approximately 40kb apart, they can be used to verify assemblies. Using Nucmer, the current reference genome was aligned to a de novo assembly constructed from Illumina and PacBio reads. Using the Nucmer alignment output, the coordinates of the fosmid alignments can be mapped from the reference to the new assembly. This allows the problem to be represented as a graph optimization problem. Each alignment block in the assembly will be represented as a node, with nodes for each contig chained together in a linked list. Edges will be constructed between/within nodes based on where each of the 2 ends of a given fosmid paired end read align. Every edge has a weight; the closer its length is to the 40kb ideal, the lower the weight. Edges between nodes on different contigs will be given an extremely high weight, since this indicates that one of the two nodes must be on the wrong contig. Using this structure, an optimization problem naturally appears, with the goal of minimizing edge weights. 

<h1>Files and Folders</h1>

block\_node.py- the class for the Node object, to support the graph data structure

delta\_parser.py- Nucmer alignment output is not is a format that can be easily used in my program. Example output can be seen [here](http://mummer.sourceforge.net/manual/#nucmer) (scroll down to the section titled The "delta" file format). This script is meant to parse this output into a form another script can use to map the coordinates of fosmid ends from the reference to the assembly.

fosmid\_edge.py- the class for the Edge object, to support the graph data structure. Note that this class currently contains the logic that determines the weight of a particular edge.

node\_list\_generator.py- main script, under development. Currently constructs all nodes and edges; however, mapping from reference to assembly coordinates has not yet been implemented, so edge coordinates are incorrect.
