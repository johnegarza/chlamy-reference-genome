import sys, os, csv, re
from block_node import Node
from fosmid_edge import Edge
import argparse
from contig_loc import ContigLocation
import pysam

############# TO RUN ###############
#	from graph directory:
#	python main.py full_node_list.txt full_unique_mapped_ends.txt

parser = argparse.ArgumentParser(description='Create initial graph from list of alignments and fosmid paired ends, then move nodes with edges in different scaffolds')
parser.add_argument('alignments', metavar='*.tsv', help='a tab delimited reference to assembly alignment file from nucmer')
parser.add_argument('fosmids', metavar='*.tsv', help = 'a fosmid paired end file')
args = parser.parse_args()

alignment_file = args.alignments
fosmid_pairs = args.fosmids

scaffolds = [] #contains head node for each scaffold
line_indexed_nodes = [] #to retrieve node at line n, call line_indexed_nodes[n-1]
bad_edges = [] 	#used for iterating over bad edges and moving them

with open(alignment_file) as a_f:
	alignment_data = csv.reader(a_f, delimiter="\t")

	#TODO refactor to use Scaffold (upcoming) data structure

	#initialize some dummy data so my script can dynamically begin construction from the 
	#beginning of the file, instead of needing to hard code the initial case
	prev_scaf = " "
	dummy_CL = ContigLocation("dummy", 0, 0)
	dummy_node = Node(-1, dummy_CL, dummy_CL)
	curr_node = dummy_node

	for line_id, block in enumerate(alignment_data):

		print(line_id)

		#load in data from the current line
		asm_scaf = block[0]
		asm_start = int(block[1])
		asm_stop = int(block[2])
		ref_chr = block[3]
		ref_start = int(block[4])
		ref_stop = int(block[5])
		line_num = int(block[8])

		if(prev_scaf != asm_scaf): #end of an assembly scaffold

			curr_node.next = None #TODO may not be necessary

			curr_ref_CL = ContigLocation(ref_chr, ref_start, ref_stop)
			curr_asm_CL = ContigLocation(asm_scaf, asm_start, asm_stop)
			curr_node = Node(line_num, curr_ref_CL, curr_asm_CL)

			#TODO again, hardcoded while testing; make dynamic parameter in later version
			region = asm_scaf + ":" + str(asm_start) + "-" + str(asm_stop)
			#faidx returns a string; it begins with the scaffold name, followed by a \n, after which one or more lines
			#of DNA are present, also separated by \n; only care about raw string sequence, hence split and [1:]
			node_seq = "".join(pysam.faidx("assembly.fasta", region).split()[1:])
			curr_node.seq = node_seq

			assert len(node_seq) == len(curr_node.asm)

			#within this case, curr_node is the first node in a scaffold, so add it to the list of head nodes
			scaffolds.append(curr_node)

			line_indexed_nodes.append(curr_node)

			prev_scaf = asm_scaf
			
			#TODO isn't this just an else case ~~~ 6/7/18- why not?
			continue #work is done for this cycle, and thanks to python's scope (or lack thereof)
				 #we can just skip to the next iteration of the loop

		new_ref_CL = ContigLocation(ref_chr, ref_start, ref_stop)
		new_asm_CL = ContigLocation(asm_scaf, asm_start, asm_stop)
		new_node = Node(line_num, new_ref_CL, new_asm_CL, p_node = curr_node)

		#TODO again, hardcoded while testing; make dynamic parameter in later version
		region = asm_scaf + ":" + str(asm_start) + "-" + str(asm_stop)
		#this returns scaffold name and sequence string separated by \n; only care about seq, hence split and [1]
		node_seq = "".join(pysam.faidx("assembly.fasta", region).split()[1:])
		new_node.seq = node_seq

		assert len(node_seq) == len(new_node.asm)

		curr_node.next = new_node
		curr_node = new_node

		line_indexed_nodes.append(curr_node)

		prev_scaf = asm_scaf

with open(fosmid_pairs) as f_p:

	pair_data = csv.reader(f_p, delimiter="\t")

	for pair in pair_data:


		left_ref_start = int(pair[1])
		left_ref_stop = int(pair[2])
		left_asm_start = int(pair[4])
		left_asm_stop = int(pair[5])
		left_block_line = int(pair[7])

		right_ref_start = int(pair[9])
		right_ref_stop = int(pair[10])
		right_asm_start = int(pair[12])
		right_asm_stop = int(pair[13])
		right_block_line = int(pair[15])


		node1 = line_indexed_nodes[left_block_line - 1]
		node2 = line_indexed_nodes[right_block_line - 1]

		edge = Edge(node1, node2, left_ref_start, left_ref_stop, left_asm_start, left_asm_stop, right_ref_start, right_ref_stop, right_asm_start, right_asm_stop)

		if node1 is node2: #prevent duplicate edges in the same node
			node1.add_edge(edge)
		else:
			node1.add_edge(edge)
			node2.add_edge(edge)

		if edge.weight == -10:
			bad_edges.append(edge)

#don't need this anymore- clear memory and unnecessary references that may keep nodes removed from main assembly alive and "orphaned"
line_indexed_nodes = []

#TODO hardcoded during testing; make this an input parameter so it's dynamic TODO
samfile = pysam.AlignmentFile("../novoalign/imp3.merged.sorted.bam", "rb")


while bad_edges: #run as long as bad_edges is not empty

	print(len(bad_edges))

	seed_edge = bad_edges[0]

	bad_node = seed_edge.node1
	other_node = seed_edge.node2


	continue_search = True
	curr_node = bad_node
	searched_nodes = [bad_node]
	region_lo = seed_edge.edge_low()
	region_hi = seed_edge.edge_high()
	left_region_node = None
	right_region_node = None

	#hack to make the loop handle all the work, even the initial case, despite possible edge case (if bad_node has sorted edges at either end
	#that match to another scaffold's node, aside from other_node, search would fail to extend boundaries from seed_edge)
	initial_edges = bad_node.get_sorted_edges()
	seed_pos = initial_edges.index(seed_edge)
	initial_left_list = initial_edges[:seed_pos]
	seed_pos += 1
	initial_right_list = initial_edges[seed_pos:]
	left_first = True
	right_first = True

	#search right
	while curr_node is not None and continue_search:

		if right_first:
			search_space = initial_right_list
			right_first = False
		else:
			search_space = curr_node.get_sorted_edges()

		for edge in search_space:

			#changed to the existing implementation because of a small but possibly significant edge case:
			#edges are sorted in ascending order of .edge_low(), which makes no guarantees as to the ordering of .edge_high()
			#for example, the following list is a valid example (with edges displayed as a tuple of (edge_low, edge_high)
			# [(1,10), (2,8), (3, 20)]
			# region_hi is intended to capture 20, but the line below would only capture up to 10 then halt
#			if (edge.weight != -10 or (edge.opposite_node(bad_node) is other_node) ) and (edge.edge_high() > region_hi):

			if edge.weight != -10 or (edge.opposite_node(bad_node) is other_node):
				if edge.edge_high() > region_hi:
					region_hi = edge.edge_high()
					right_region_node = curr_node
			else:
				#break out of this while loop
				continue_search = False
				break


		#note that these 2 lines will always run- so the first and last nodes in searched_nodes will be one prior and one after
		#the last node searched at that extreme; this could potentially be NoneType, allowing for easy identification of scaffold-spanning
		#chunk regions later
		curr_node = curr_node.next
		searched_nodes.append(curr_node)



	continue_search = True
	curr_node = bad_node
	#want to be able to efficiently add to the beginning of this list; still want efficient random access, so no deque (TODO- is random access necessary?)
	searched_nodes.reverse()

	while curr_node is not None and continue_search:

		if left_first:
			search_space = initial_left_list
			left_first = False
		else
			search_space = curr_node.get_sorted_edges()

		search_space.reverse()

		for edge in search_space:
			if (edge.weight != -10 or (edge.opposite_node(bad_node) is other_node) ) and (edge.edge_low() < region_lo):
				region_lo = edge.edge_low()
				left_region_node = curr_node
			else:
				continue_search = False
				break

		curr_node = curr_node.prev
		searched_nodes.append(curr_node)

	searched_nodes.reverse() #back to normal

	print(len(searched_nodes))

