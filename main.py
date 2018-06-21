from collections import defaultdict
import sys, os, csv, re, fnmatch
from block_node import Node
from fosmid_edge import Edge
import pickle
import argparse
import copy
from contig_loc import ContigLocation

if len(sys.argv) < 3:
	sys.exit("Usage: %s block_list_tab_delimited indexed_fosmid_pairs" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])
if not os.path.exists(sys.argv[2]):
	sys.exit("Error: File '%s' not found" % sys.argv[2])

############# TO RUN ###############
#	from graph directory:
#	python main.py htcf_data/node_list.tsv unique_mapped_ends.txt

alignment_file = sys.argv[1] #tab_delim_results.tsv
fosmid_pairs = sys.argv[2]

'''
parser = argparse.ArgumentParser(description='creates nodes and edges, then associates them together into a graph')
parser.add_argument('alignments', metavar='*.tsv', help='a tab delimited reference to assembly alignment file')
parser.add_argument('fosmids', metavar='*.tsv', help = 'a fosmid paired end file')
args = parser.parse_args()

alignment_file = args.alignments #tab_delim_results.tsv
fosmid_pairs = args.fosmids
'''

contigs = [] #contains head node for each contig
line_indexed_nodes = [] #to retrieve node at line n, call line_indexed_nodes[n-1]
#edges = [] #adding to help with active development
bad_edges = [] 	#will attempt to use this for iterating over bad edges and moving them, instead of 
		#iterating over nodes and searching each one's list of edges

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

		#load in data from the current line
		asm_scaf = block[0]
		asm_start = int(block[1])
		asm_stop = int(block[2])
		ref_chr = block[3]
		ref_start = int(block[4])
		ref_stop = int(block[5])
		line_num = int(block[8])

		if(prev_scaf != asm_scaf): #end of an assembly contig

			curr_node.next = None #TODO may not be necessary

			curr_ref_CL = ContigLocation(ref_chr, ref_start, ref_stop)
			curr_asm_CL = ContigLocation(asm_scaf, asm_start, asm_stop)
			curr_node = Node(line_num, curr_ref_CL, curr_asm_CL)

			contigs.append(curr_node)

			line_indexed_nodes.append(curr_node)

			prev_scaf = asm_scaf
			
			#TODO isn't this just an else case ~~~ 6/7/18- why not?
			continue #work is done for this cycle, and thanks to python's scope (or lack thereof)
				 #we can just skip to the next iteration of the loop

		new_ref_CL = ContigLocation(ref_chr, ref_start, ref_stop)
		new_asm_CL = ContigLocation(asm_scaf, asm_start, asm_stop)
		new_node = Node(line_num, new_ref_CL, new_asm_CL, p_node = curr_node)

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
		#edges.append(edge) #TODO if no proper use for this, remove; will just lead to memory leaks, as this keeps edges deleted later on still alive due to the reference

'''
tally = [0]*200
for cnum, contig_head in enumerate(contigs):
	
	iterator = contig_head
	while(iterator is not None):
		counter = 0
		for edge in iterator._edges:
			if edge.weight == -10:
				counter += 1
#		if counter != 0:
#			print(str(counter) + "/" + str(len(iterator._edges)) + " bad edges")
		tally[counter] += 1
		iterator = iterator.next

for place, val in enumerate(tally):
	if val != 0:
		print( str(place) + " edges: " + str(val) + " nodes" )
'''

#don't need this anymore- clear memory and unnecessary references that may keep nodes removed from main assembly alive and "orphaned"
lined_indexed_nodes = []

debug_index = 0
while bad_edges: #run as long as bad_edges is not empty

	seed_edge = bad_edges[0]

	#arbitrarily chose to start with node1; will work with node2 later
	bad_node = seed_edge.node1
	other_node = seed_edge.node2
	assert(bad_node is not other_node) #by definition this should be correct

	### STEP 1 Group up edges ###

	search_space = bad_node.get_sorted_edges()
	assert(len(search_space) > 0)

	if seed_edge not in search_space:
		print(str(seed_edge))
		print(str(seed_edge.node1))
		print(id(seed_edge.node1))
		print(str(seed_edge.node2))
		print(id(seed_edge.node2))
		print(str(bad_node))
		print(id(bad_node))
	assert(seed_edge in search_space)


	print("round " + str(debug_index) + " / " + str(len(bad_edges)) )
#	for edge in search_space:
#		print(str(edge))

	chunk_lo = float('inf')
	chunk_hi = float('-inf')
	
	for edge in search_space:

		if edge.edge_low(bad_node) < chunk_lo and edge.opposite_node(bad_node) is other_node:
			chunk_lo = edge.edge_low(bad_node)

		if edge.edge_high(bad_node) > chunk_hi and edge.opposite_node(bad_node) is other_node:
			chunk_hi = edge.edge_high(bad_node)

	left_edges = []
	left_border_edges = []
	chunk_edges = []
	right_border_edges = []
	right_edges = []

	stop = len(search_space)
	index = 0
	curr_edge = search_space[0]

	#TODO note: in the following 2 while loops, index < stop and index >= stop may be redundant, if needed at all
	#second condition shouldn't be necessary; shouldn't affect performance, but may reevaluate later
	while ( curr_edge.edge_low(bad_node) < chunk_lo ) and (index < stop):

		if curr_edge.edge_high(bad_node) <= chunk_lo:
			left_edges.append(curr_edge)
		else:
			left_border_edges.append(curr_edge)

		index += 1
		if index >= stop: #for an edge case that may not be possible; reevaluate later
			break
		curr_edge = search_space[index]

	#second condition shouldn't be necessary; shouldn't affect performance, but may reevaluate later
	while ( curr_edge.edge_low(bad_node) < chunk_hi ) and (index < stop):

		if curr_edge.edge_high(bad_node) <= chunk_hi:
			chunk_edges.append(curr_edge)
		else:
			right_border_edges.append(curr_edge)

		index += 1
		if index >= stop: #for an edge case that may not be possible; reevaluate later
			break
		curr_edge = search_space[index]

	while index < stop:
		right_edges.append(curr_edge)
		index += 1
		if index >= stop:
			break
		curr_edge = search_space[index]

	total_len = len(left_edges) + len(left_border_edges) + len(chunk_edges) + len(right_border_edges) + len(right_edges)
	assert( total_len == stop )
#	print("good")
#	print( str(len(chunk_edges)) + " " + str(stop) )

	if seed_edge not in chunk_edges:
		print(str(seed_edge))
		print(str(chunk_lo) + "-" + str(chunk_hi))
	assert(seed_edge in chunk_edges)

	assert(seed_edge not in left_border_edges)
	assert(seed_edge not in right_border_edges)

	b_e_temp_set = set(bad_edges)
	b_e_temp_set.difference_update(chunk_edges) #anything in chunk_edges that's also in b_e will be removed from b_e
	bad_edges = list(b_e_temp_set)

	print(len(chunk_edges))

	'''
	if (len(left_border_edges) != 0):
		print("lb edges: " + str(len(left_border_edges)))

	if (len(right_border_edges) != 0):
		print("rb edges: " + str(len(right_border_edges)))
	'''
	
	#TODO these edges will be important in the real algorithm, but for now while testing basic ops these will just be ignored
	b_e_temp_set = set(bad_edges)
	b_e_temp_set.difference_update(left_border_edges)
	b_e_temp_set.difference_update(right_border_edges)
	bad_edges = list(b_e_temp_set)

	for edge in left_border_edges:
		bad_node.remove_edge(edge)
		temp = edge.opposite_node(bad_node)
		if temp is not bad_node: #fix removal failures when both edge endpoints are in the same node
			temp.remove_edge(edge)
	left_border_edges = []

	for edge in right_border_edges:
		bad_node.remove_edge(edge)
		temp = edge.opposite_node(bad_node)
		if temp is not bad_node: #fix removal failures when both edge endpoints are in the same node
			temp.remove_edge(edge)
	right_border_edges = []
	#end TODO	

	######################## CREATE NEW NODES ####################
	'''
	node_len = (chunk_hi - chunk_lo) + 1

	left_dist = chunk_lo - bad_node.asm.left #inclusive coords
	right_dist = bad_node.asm.right - chunk_hi


	#TODO ref CL changes only work if nodes are a one-to-one mapping; is this accurate?

	chunk_ref_CL = bad_node.ref.trim(left_dist, right_dist)
	chunk_asm_CL = ContigLocation(other_node.asm.name, other_node.asm.left, other_node.asm.left + (node_len - 1) )
	chunk_node = Node(-1, chunk_ref_CL, chunk_asm_CL, bad_node.asm_original, chunk_edges)

	#TODO should this be      .trim_right?
	left_ref_CL = bad_node.ref.trim_left(left_dist - 1) #-1 because otherwise this and prev node would start at the exact same coord; this CL should have exclusive coords
	left_asm_CL = ContigLocation(bad_node.asm.name, bad_node.asm.left, chunk_lo - 1)
	left_node = Node(-1, left_ref_CL, left_asm_CL, bad_node.asm_original, left_edges)

	#TODO should this be       .trim_left?
	right_ref_CL = bad_node.ref.trim_right(right_dist - 1)
	right_asm_CL = ContigLocation(bad_node.asm.name, chunk_lo, bad_node.asm.right - node_len)
	right_node = Node(-1, right_ref_CL, right_asm_CL, bad_node.asm_original, right_edges)

	full_len = bad_node.asm.right - bad_node.asm.left
	nodes_len = (left_node.asm.right - left_node.asm.left) + (chunk_node.asm.right - chunk_node.asm.left) + (right_node.asm.right - right_node.asm.left) + 2
	assert(full_len == nodes_len)
	'''
	assert(seed_edge in bad_node.get_sorted_edges())

	debug_index += 1



