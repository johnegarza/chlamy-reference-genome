from collections import defaultdict
import sys, os, csv, re, fnmatch
from block_node import Node
from fosmid_edge import Edge
import pickle
import argparse
import copy
from contig_loc import ContigLocation
import gc
import pysam

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

		print(line_id)

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

			#TODO again, hardcoded while testing; make dynamic parameter in later version
			region = asm_scaf + ":" + str(asm_start) + "-" + str(asm_stop)
			#this returns scaffold name and sequence string separated by \n; only care about seq, hence split and [1:]
			node_seq = "".join(pysam.faidx("assembly.fasta", region).split()[1:])
			curr_node.seq = node_seq

			assert len(node_seq) == len(curr_node.asm)

			contigs.append(curr_node)

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
#			assert(node1.asm.left <= left_asm_start)
#			assert(node1.asm.right >= left_asm_stop)
#			assert(node1.asm.left <= right_asm_start)
#			assert(node1.asm.right >= right_asm_stop)
		else:
			node1.add_edge(edge)
			node2.add_edge(edge)
#			assert(node1.asm.left <= left_asm_start)
#			assert(node1.asm.right >= left_asm_stop)
#			assert(node2.asm.left <= right_asm_start)
#			assert(node2.asm.right >= right_asm_stop)


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
line_indexed_nodes = []

#TODO hardcoded during testing; make this an input parameter so it's dynamic TODO
samfile = pysam.AlignmentFile("../novoalign/imp3.merged.sorted.bam", "rb")


while bad_edges: #run as long as bad_edges is not empty

	print(len(bad_edges))

	seed_edge = bad_edges[0]

	#arbitrarily chose to start with node1; will work with node2 later
	bad_node = seed_edge.node1
	other_node = seed_edge.node2

	### STEP 1 Group up edges ###

	search_space1 = bad_node.get_sorted_edges()
	assert(search_space1[0].edge_low(bad_node) > bad_node.asm.left)
	assert(search_space1[-1].edge_high(bad_node) < bad_node.asm.right)

	chunk_lo1 = float('inf')
	chunk_hi1 = float('-inf')
	lo_pos1 = -1
	hi_pos1 = -1
	
	for pos, edge in enumerate(search_space1):

		if edge.edge_low(bad_node) < chunk_lo1 and edge.opposite_node(bad_node) is other_node:
			chunk_lo1 = edge.edge_low(bad_node)
			depth_lo1 = edge.og_edge_low(bad_node)
			lo_pos1 = pos

		if edge.edge_high(bad_node) > chunk_hi1 and edge.opposite_node(bad_node) is other_node:
			chunk_hi1 = edge.edge_high(bad_node)
			depth_hi1 = edge.og_edge_high(bad_node)
			hi_pos1 = pos

	assert(lo_pos1 != -1)
	assert(hi_pos1 != -1)
	assert(chunk_lo1 >= bad_node.asm.left)
	assert(chunk_hi1 <= bad_node.asm.right)

	#TODO this search is designed to ignore border edges; should update later to account for these (if any are present)
	lo_search1 = depth_lo1
	left_search = list(reversed(search_space1[:lo_pos1]))
	for edge in left_search:
		if edge.edge_high(bad_node) < chunk_lo1:
			lo_search1 = edge.og_edge_high(bad_node)
			break

	hi_search1 = depth_hi1
	right_search = []
	hi_pos1 += 1 #since slice will use this inclusively, and I don't want an already examined node to be included
	if hi_pos1 < len(search_space1):
		right_search = search_space1[hi_pos1:]
	for edge in right_search:
		if edge.edge_low(bad_node) > chunk_hi1:
			hi_search1 = edge.og_edge_low(bad_node)
			break
	
	lo_search1 += 1 #to make these inclusive coords for the depth of coverage search space
	hi_search1 -= 1

	#########################################################################################################################################
	search_space2 = other_node.get_sorted_edges()
	assert(search_space2[0].edge_low(other_node) > other_node.asm.left)
	assert(search_space2[-1].edge_high(other_node) < other_node.asm.right)

	chunk_lo2 = float('inf')
	chunk_hi2 = float('-inf')
	lo_pos2 = -1
	hi_pos2 = -1
	
	for pos, edge in enumerate(search_space2):

		if edge.edge_low(other_node) < chunk_lo2 and edge.opposite_node(other_node) is bad_node:
			chunk_lo2 = edge.edge_low(other_node)
			depth_lo2 = edge.og_edge_low(other_node)
			lo_pos2 = pos

		if edge.edge_high(other_node) > chunk_hi2 and edge.opposite_node(other_node) is bad_node:
			chunk_hi2 = edge.edge_high(other_node)
			depth_hi2 = edge.og_edge_high(other_node)
			hi_pos2 = pos

	assert(lo_pos2 != -1)
	assert(hi_pos2 != -1)
	assert(chunk_lo2 >= other_node.asm.left)
	assert(chunk_hi2 <= other_node.asm.right)

	#TODO this search is designed to ignore border edges; should update later to account for these (if any are present)
	lo_search2 = depth_lo2
	left_search = list(reversed(search_space2[:lo_pos2]))
	for edge in left_search:
		if edge.edge_high(other_node) < chunk_lo2:
			lo_search2 = edge.og_edge_high(other_node)
			break

	hi_search2 = depth_hi2
	right_search = []
	hi_pos2 += 1 #since slice will use this inclusively, and I don't want an already examined node to be included
	if hi_pos2 < len(search_space2):
		right_search = search_space2[hi_pos2:]
	for edge in right_search:
		if edge.edge_low(other_node) > chunk_hi2:
			hi_search2 = edge.og_edge_low(other_node)
			break
	
	lo_search2 += 1 #to make these inclusive coords for the depth of coverage search space
	hi_search2 -= 1

	###############################################
	#choose new chunk_lo, chunk_hi, and potentially swap bad_node/other_node based on the search region with the lowest depth of coverage

	left1_check_region = []
	right1_check_region = []
	left2_check_region = []
	right2_check_region = []

	if ( (depth_lo1-1) - lo_search1 ) >= 0:
		temp = samfile.pileup(bad_node.asm_original.name, lo_search1, (depth_lo1 - 1) )
		left1_check_region = [ x for x in temp if lo_search1 <= x.pos < depth_lo1 ]

	if ( hi_search1 - (depth_hi1+1) ) >= 0:
		temp = samfile.pileup(bad_node.asm_original.name, (depth_hi1+1), hi_search1 )
		right1_check_region = [ x for x in temp if depth_hi1 < x.pos <= hi_search1 ]

	if ( (depth_lo2-1) - lo_search2 ) >= 0:
		temp = samfile.pileup(other_node.asm_original.name, lo_search2, (depth_lo2 - 1) )
		left2_check_region = [ x for x in temp if lo_search2 <= x.pos < depth_lo2 ]

	if ( hi_search2 - (depth_hi2+1) ) >= 0:
		temp = samfile.pileup(other_node.asm_original.name, (depth_hi2+1), hi_search2 )
		right2_check_region = [ x for x in temp if depth_hi2 < x.pos <= hi_search2 ]

	assert(left2_check_region is not None)

	assert(left2_check_region is not None)

	region1_len = len(left1_check_region) + len(right1_check_region)
	if region1_len != 0:
		avg1 = float( sum( [x.pos for x in left1_check_region] ) + sum( [x.pos for x in right1_check_region] ) ) / region1_len
	else:
		avg1 = float('inf')
	region2_len = len(left2_check_region) + len(right2_check_region)

	if region2_len != 0:
		avg2 = float( sum( [x.pos for x in left2_check_region] ) + sum( [x.pos for x in right2_check_region] ) ) / region2_len
	else:
		avg2 = float('inf')

	if avg1 < avg2: #will move bad_node to other_node's scaffold

		search_space = search_space1 #set parameter for rest of move algorithm

		min_pos = ( float('inf'), None )
		left1_check_region.reverse() #reverse so we're conservative and find the minimum closer to chunk_lo1 instead of lo_search1
		for pos in left1_check_region:
			if pos.n < min_pos[0]:
				min_pos = (pos.n, pos.pos)

		if min_pos[1] is not None: #set parameter for rest of move algorithm
			left_diff = depth_lo1 - min_pos[1]
			chunk_lo = chunk_lo1 - left_diff 
		else:
			chunk_lo = chunk_lo1


		min_pos = ( float('inf'), None )
		for pos in right1_check_region:
			if pos.n < min_pos[0]:
				min_pos = (pos.n, pos.pos)

		if min_pos[1] is not None: #set parameter for rest of move algorithm
			right_diff = min_pos[1] - depth_hi1
			chunk_hi = chunk_hi1 + right_diff 
		else:
			chunk_hi = chunk_hi1



	else: #will move other_node to bad_node's scaffold
		search_space = search_space2

		min_pos = ( float('inf'), None )
		left2_check_region.reverse() #reverse so we're conservative and find the minimum closer to chunk_lo1 instead of lo_search1
		for pos in left2_check_region:
			if pos.n < min_pos[0]:
				min_pos = (pos.n, pos.pos)

		if min_pos[1] is not None: #set parameter for rest of move algorithm
			left_diff = depth_lo2 - min_pos[1]
			chunk_lo = chunk_lo2 - left_diff 
		else:
			chunk_lo = chunk_lo2


		min_pos = ( float('inf'), None )
		for pos in right2_check_region:
			if pos.n < min_pos[0]:
				min_pos = (pos.n, pos.pos)

		if min_pos[1] is not None: #set parameter for rest of move algorithm
			right_diff = min_pos[1] - depth_hi2
			chunk_hi = chunk_hi2 + right_diff 
		else:
			chunk_hi = chunk_hi2

		temp = bad_node
		bad_node = other_node
		other_node = temp

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

	for edge in left_edges:
		assert(edge.edge_low(bad_node) > bad_node.asm.left)
		assert(edge.edge_high(bad_node) < bad_node.asm.right)
	for edge in left_border_edges:
		assert(edge.edge_low(bad_node) > bad_node.asm.left)
		assert(edge.edge_high(bad_node) < bad_node.asm.right)

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

	for edge in chunk_edges:
		assert(edge.edge_low(bad_node) > bad_node.asm.left)
		assert(edge.edge_high(bad_node) < bad_node.asm.right)

	for edge in right_border_edges:
		assert(edge.edge_low(bad_node) > bad_node.asm.left)
		assert(edge.edge_high(bad_node) < bad_node.asm.right)


	while index < stop:
		right_edges.append(curr_edge)
		index += 1
		if index >= stop:
			break
		curr_edge = search_space[index]

	for edge in right_edges:
		assert(edge.edge_low(bad_node) > bad_node.asm.left)
		assert(edge.edge_high(bad_node) < bad_node.asm.right)

	total_len = len(left_edges) + len(left_border_edges) + len(chunk_edges) + len(right_border_edges) + len(right_edges)
	assert( total_len == stop )

	b_e_temp_set = set(bad_edges)
	b_e_temp_set.difference_update(chunk_edges) #anything in chunk_edges that's also in b_e will be removed from b_e
	bad_edges = list(b_e_temp_set)

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
	
	node_len = (chunk_hi - chunk_lo) + 1

	left_dist = chunk_lo - bad_node.asm.left #inclusive coords
	right_dist = bad_node.asm.right - chunk_hi
	right_split_index = (chunk_hi - bad_node.asm.left) + 1


	#TODO ref CL changes only work if nodes are a one-to-one mapping; is this accurate?

	chunk_ref_CL = bad_node.ref.trim(left_dist, right_dist)
	chunk_asm_CL = ContigLocation(other_node.asm.name, other_node.asm.left, other_node.asm.left + (node_len - 1) )
	chunk_node = Node(-1, chunk_ref_CL, chunk_asm_CL, bad_node.asm_original, chunk_edges)
	chunk_seq = bad_node.seq[left_dist:right_split_index]
	chunk_node.seq = chunk_seq

	assert len(chunk_seq) == len(chunk_asm_CL)

	#TODO should this be      .trim_right?
	left_ref_CL = bad_node.ref.trim_right(right_dist - 1) #-1 because otherwise this and prev node would start at the exact same coord; this CL should have exclusive coords
	left_asm_CL = ContigLocation(bad_node.asm.name, bad_node.asm.left, chunk_lo - 1)
	left_node = Node(-1, left_ref_CL, left_asm_CL, bad_node.asm_original, left_edges)
	left_seq = bad_node.seq[:left_dist]
	left_node.seq = left_seq

	assert len(left_seq) == len(left_asm_CL)

	#TODO should this be       .trim_left?
	right_ref_CL = bad_node.ref.trim_left(left_dist - 1)
	right_asm_CL = ContigLocation(bad_node.asm.name, chunk_lo, bad_node.asm.right - node_len)
	right_node = Node(-1, right_ref_CL, right_asm_CL, bad_node.asm_original, right_edges)
	right_seq = bad_node.seq[right_split_index:]
	right_node.seq = right_seq

	assert len(right_seq) == len(right_asm_CL)

	full_len = bad_node.asm.right - bad_node.asm.left
	nodes_len = (left_node.asm.right - left_node.asm.left) + (chunk_node.asm.right - chunk_node.asm.left) + (right_node.asm.right - right_node.asm.left) + 2
	assert(full_len == nodes_len)
	
	####INSERT NEW NODES####
	
	left_node.next = right_node
	right_node.prev = left_node

	#node could be the head of a contig
	if bad_node.prev is not None:
		bad_node.prev.next = left_node
		left_node.prev = bad_node.prev
	else:
		for index, head in enumerate(contigs):
			if head is bad_node:
				contigs[index] = left_node

	#node could be the tail of a contig
	if bad_node.next is not None:
		bad_node.next.prev = right_node
		right_node.next = bad_node.next



	#other_node could be the head of a contig
	if other_node.prev is None:
		for index, head in enumerate(contigs):
			if head is other_node:
				contigs[index] = chunk_node
	else:
		chunk_node.prev = other_node.prev
		other_node.prev.next = chunk_node

	chunk_node.next = other_node
	other_node.prev = chunk_node

	###### UPDATE EDGE ENDPOINTS ########
	left_node.new_edge_endpoints(bad_node)
	right_node.new_edge_endpoints(bad_node)
	chunk_node.new_edge_endpoints(bad_node)
	bad_node.clear()

	
	#STEP 6 propogate coordinate updates to all nodes (and their edges) that changed location due to the insertion/deletion

	#TODO refactor so that chunk_node is part of the loop? removes many of the next lines and makes this cleaner overall
	chunk_node.shift_edges(node_len) #coords have already been shifted during creation
	if chunk_node.next is not None:
		iterator = chunk_node.next
		if iterator.next is None:
			iterator.shift(node_len)
		while iterator.next is not None:
			iterator.shift(node_len)
			iterator = iterator.next
		
	#TODO same as above
	right_node.shift_edges(-node_len)
	if right_node.next is not None:
		iterator = right_node.next
		if iterator.next is None:
			iterator.shift(-node_len)
		while iterator.next is not None:
			iterator.shift(-node_len)
			iterator = iterator.next

	for edge in left_edges:
		assert(edge.edge_low(left_node) > left_node.asm.left)
		assert(edge.edge_high(left_node) < left_node.asm.right)
	for edge in chunk_edges:
		assert(edge.edge_low(chunk_node) > chunk_node.asm.left)
		assert(edge.edge_high(chunk_node) < chunk_node.asm.right)
	for edge in right_edges:
		assert(edge.edge_low(right_node) > right_node.asm.left)
		assert(edge.edge_high(right_node) < right_node.asm.right)


samfile.close()

'''	
for head in contigs:
	while head is not None:
		print(head)
		head = head.next
	print("\n")
'''

with open("output2.txt", "w") as o_f:
	for head in contigs:
		while head is not None:
			o_f.write(str(head))
			o_f.write("\n")
			head = head.next

		o_f.write("\n")







