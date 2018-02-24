#The Problem#

Experimental evidence has demonstrated that at least one large portion of the C. reinhardtii genome is misaligned in the current reference genome, with others suspected. The problem has persisted across multiple drafts of the genome.

#The Solution (?)#
A large library of BAC and fosmid paired-end reads aligned to the reference genome is available on phytozome. This dataset was downloaded and processed to remove the BAC data. Since a very high percentage of fosmid paired ends are approximately 40kb apart, they can be used to verify assemblies. Using Nucmer, the current reference genome was aligned to a de novo assembly constructed from Illumina and PacBio reads. Using the Nucmer alignment output, the coordinates of the fosmid alignments can be mapped from the reference to the new assembly. This allows the problem to be represented as a graph optimization problem. Each alignment block in the assembly will be represented as a node, with nodes for each contig chained together in a linked list. Edges will be constructed between/within nodes based on where each of the 2 ends of a given fosmid paired end read align. Every edge has a weight; the closer its length is to the 40kb ideal, the lower the weight. Edges between nodes on different contigs will be given an extremely high weight, since this indicates that one of the two nodes must be on the wrong contig. Using this structure, an optimization problem naturally appears, with the goal of minimizing edge weights. 

#Files and Folders#
