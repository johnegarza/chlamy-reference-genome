import sys, os, csv, re
from block_node import Node
from fosmid_edge import Edge
import argparse
from contig_loc import ContigLocation
import pysam
import pdb

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

debug_trimmed_nodes = set()

debug_bad_edge = None
debug_edge_info = []
while bad_edges: #run as long as bad_edges is not empty

	print(len(bad_edges))


	if len(bad_edges) == 1161:
		#offending edge is at position bad_edges[1]

		debug_bad_edge = bad_edges[1]
		debug_edge_info.append("set d_b_e")

		'''
		spots = []
		for index, edge in enumerate(bad_edges):
			if edge.asm1_original.name == "chromosome_3" and edge.asm1_original.left == 187189:
				spots.append(index)
		import pdb; pdb.set_trace()
		'''

	seed_edge = bad_edges[0]

	bad_node = seed_edge.node1
	other_node = seed_edge.node2
	insert_left = other_node.prev

	assert seed_edge.edge_low(bad_node) >= bad_node.asm.low()
	assert seed_edge.edge_high(bad_node) <= bad_node.asm.high()
	assert seed_edge.edge_low(other_node) >= other_node.asm.low()
	assert seed_edge.edge_high(other_node) <= other_node.asm.high()

	'''
	######## define the low/high coordinates/nodes of the region to be pulled out of its current scaffold and placed in another ################
	'''

	continue_search = True
	curr_node = bad_node
	searched_nodes = [bad_node]
	region_lo = seed_edge.edge_low(bad_node)
	region_hi = seed_edge.edge_high(bad_node)
	left_region_node = None
	right_region_node = None

	#hack to make the loop handle all the work, even the initial case, despite possible edge case (if bad_node has sorted edges at either end
	#that match to another scaffold's node, aside from other_node, search would fail to extend boundaries from seed_edge)
	initial_edges = bad_node.get_sorted_edges()

	for edge in initial_edges:
		if edge is debug_bad_edge:
			debug_edge_info.append(189)
			
		assert edge.edge_low(bad_node) >= bad_node.asm.low()
		assert edge.edge_high(bad_node) <= bad_node.asm.high()
	'''
	seed_pos = initial_edges.index(seed_edge)
	initial_left_list = initial_edges[:seed_pos]
	seed_pos += 1
	initial_right_list = initial_edges[(seed_pos + 1):]
	'''

	initial_left_list = []
	initial_right_list = []
	for edge in initial_edges:
		if edge.edge_high(bad_node) <= seed_edge.edge_low(bad_node):
			initial_left_list.append(edge)
		elif edge.edge_low(bad_node) >= seed_edge.edge_high(bad_node):
			initial_right_list.append(edge)

	left_first = True
	right_first = True

	edge_list = set()
	edge_list.add(seed_edge)

	#search right
	right_debug = bad_node
	right_first_iter = True
	while curr_node is not None and continue_search:
		assert(curr_node is not other_node)
		if right_first:
			search_space = initial_right_list
			for edge in search_space:
				assert edge.edge_low(curr_node) >= seed_edge.edge_high(curr_node)
				assert edge.edge_high(curr_node) <= curr_node.asm.high()
			right_first = False
		else:
			search_space = curr_node.get_sorted_edges()

		right_exclusive_edges = []

		for edge in search_space:
			if edge is debug_bad_edge:
				debug_edge_info.append((232, "right scan", curr_node))
			assert edge.edge_low(curr_node) >= curr_node.asm.low()
			assert edge.edge_high(curr_node) <= curr_node.asm.high()
			assert edge.node1 is curr_node or edge.node2 is curr_node


		node_searched = False
		for edge in search_space:

			#changed to the existing implementation because of a small but possibly significant edge case:
			#edges are sorted in ascending order of .edge_low(), which makes no guarantees as to the ordering of .edge_high()
			#for example, the following list is a valid example (with edges displayed as a tuple of (edge_low, edge_high)
			# [(1,10), (2,8), (3, 20)]
			# region_hi is intended to capture 20, but the line below would only capture up to 10 then halt
#			if (edge.weight != -10 or (edge.opposite_node(bad_node) is other_node) ) and (edge.edge_high() > region_hi):
#			import pdb; pdb.set_trace()
			if edge.weight != -10 or (edge.opposite_node(curr_node) is other_node):
				if edge.edge_high(curr_node) > region_hi:
					node_searched = True
					region_hi = edge.edge_high(curr_node)
					DEBUG_RIGHT_EDGE = edge
					assert curr_node.asm.low() <= region_hi <= curr_node.asm.high()
					right_region_node = curr_node
					edge_list.add(edge)
					right_exclusive_edges.append(edge)
			else:
				#break out of this while loop
				continue_search = False
				break

#		right_node_exists = curr_node.asm.high() - region_hi > 1 BROKEN FOR UNKNOWN REASONS

		if node_searched or right_first_iter:
			right_debug = curr_node
			#note that these 2 lines will always run- so the first and last nodes in searched_nodes will be one prior and one after
			#the last node searched at that extreme; this could potentially be NoneType, allowing for easy identification of scaffold-spanning
			#chunk regions later
			curr_node = curr_node.next
			searched_nodes.append(curr_node)
		else:
			break
		right_first_iter = False

	assert(other_node not in searched_nodes)

	continue_search = True
	curr_node = bad_node
	#want to be able to efficiently add to the beginning of this list; still want efficient random access, so no deque (TODO- is random access necessary?)
	searched_nodes.reverse()

	#search left
	left_debug = bad_node
	left_first_iter = True
	while curr_node is not None and continue_search:
		assert(curr_node is not other_node)
		if left_first:
			search_space = initial_left_list
			for edge in search_space:
				assert edge.edge_low(curr_node) >= curr_node.asm.low()
				assert edge.edge_high(curr_node) <= seed_edge.edge_low(curr_node)
			left_first = False
		else:
			search_space = curr_node.get_sorted_edges()

		for edge in search_space:
			if edge is debug_bad_edge:
				debug_edge_info.append((232, "right scan", curr_node))

			assert edge.edge_low(curr_node) >= curr_node.asm.low()
			assert edge.edge_high(curr_node) <= curr_node.asm.high()

		#otherwise this would always choose search_space[0] due to ordering
		search_space.reverse()

		left_exclusive_edges = []

		node_searched = False
		for edge in search_space:
			if (edge.weight != -10 or (edge.opposite_node(curr_node) is other_node) ) and (edge.edge_low(curr_node) < region_lo):
				node_searched = True
				region_lo = edge.edge_low(curr_node)
				DEBUG_LEFT_EDGE = edge
				assert curr_node.asm.low() <= region_lo <= curr_node.asm.high()
				left_region_node = curr_node
				edge_list.add(edge)
				left_exclusive_edges.append(edge)
			else:
				continue_search = False
				break

#		left_node_exists = region_lo - curr_node.asm.low() > 1 possibly also broken; refactoring both for consistency

		#if node_searched is still false, no values were updated based on this node; it is a fringe node and would have
		#already been added in the last iteration of the while loop, so nothing needs to be added to searched_nodes
		if node_searched or left_first_iter:
			left_debug = curr_node
			curr_node = curr_node.prev
			searched_nodes.append(curr_node)
		else:
			break
		left_first_iter = False

	searched_nodes.reverse() #back to normal
	assert(other_node not in searched_nodes)

#	if (len(searched_nodes)-2 > 20 ):
#		print(region_hi - region_lo)

	debug_trimmed_nodes.add(left_debug)
	debug_trimmed_nodes.add(right_debug)

	assert(left_debug is searched_nodes[1])
	assert(right_debug is searched_nodes[-2])

	left_node_exists = region_lo - left_debug.asm.low() > 1
	right_node_exists = region_hi - right_debug.asm.high() > 1

	assert( searched_nodes[1].asm.low() <= region_lo <= searched_nodes[1].asm.high() )
	assert( searched_nodes[-2].asm.low() <= region_hi <= searched_nodes[-2].asm.high() )

	'''
	################# begin construction of possible new nodes at either end of the region ###################
	'''

	edges_to_remove = []

	chunk_len = (region_hi - region_lo) + 1
#	right_split_index = region_hi - 

	if left_node_exists:

		left_edges = searched_nodes[1].get_sorted_edges()
		for edge in left_edges:
			assert edge.edge_low(searched_nodes[1]) >= searched_nodes[1].asm.low()
			assert edge.edge_high(searched_nodes[1]) <= searched_nodes[1].asm.high()

		left_node_edges = []
		border_edges = []
		left_chunk_edges = []
		for edge in left_edges:
			if edge.edge_low(searched_nodes[1]) < region_lo:
				if edge.edge_high(searched_nodes[1]) < region_lo:
					left_node_edges.append(edge)
				else:
					border_edges.append(edge)
			else:
				left_chunk_edges.append(edge)

		edges_to_remove.extend(border_edges)
		edges_to_remove.extend(left_chunk_edges)

		'''
		#left_exclusive_edges has all edges in searched_nodes[1] that DO NOT belong in the new left node, plus others;
		#difference update removes all elements that occur in its argument
		left_edges.difference_update(left_exclusive_edges)
		for edge in left_edges:
			assert edge.edge_low(searched_nodes[1]) >= searched_nodes[1].asm.low()
			assert edge.edge_high(searched_nodes[1]) <= (region_lo - 1)


		left_chunk_edges = set(searched_nodes[1].get_edges())
		left_chunk_edges.difference_update(left_edges)
		'''
		searched_nodes[1]._edges = left_chunk_edges
		searched_nodes[1]._edges_sorted = False

		debug_temp = searched_nodes[1].asm.low()

		left_trim_dist = region_lo - searched_nodes[1].asm.low()
		left_ref_CL = searched_nodes[1].ref.trim_lo(left_trim_dist) 
		left_asm_CL = searched_nodes[1].asm.trim_lo(left_trim_dist)
		assert( left_asm_CL.low() == debug_temp )
		assert( left_asm_CL.high() == region_lo - 1 )
		assert( searched_nodes[1].asm.low() == region_lo )
		left_asm_og_CL = searched_nodes[1].asm_original.trim_lo(left_trim_dist)
		left_node = Node(-1, left_ref_CL, left_asm_CL, left_asm_og_CL, left_node_edges)
		for edge in left_node_edges:
			#searched_nodes[1] because ownership of these edges has not been transferred from this to left_node yet
			assert( edge.edge_low(searched_nodes[1]) >= left_node.asm.low())
			assert( edge.edge_high(searched_nodes[1]) <= left_node.asm.high())
		left_seq = searched_nodes[1].seq[:left_trim_dist]
		left_node.seq = left_seq
		assert len(left_seq) == len(left_asm_CL)
		searched_nodes[1].seq = searched_nodes[1].seq[left_trim_dist:]
		assert len(searched_nodes[1].seq) == len(searched_nodes[1].asm)

	if right_node_exists:

		right_edges = searched_nodes[-2].get_sorted_edges()
#		right_edges.difference_update(right_exclusive_edges)

		right_node_edges = []
		border_edges2 = []
		right_chunk_edges = []
		for edge in right_edges:
			if edge.edge_low(searched_nodes[-2]) <= region_hi:
				if edge.edge_high(searched_nodes[-2]) <= region_hi:
					right_chunk_edges.append(edge)
				else:
					border_edges2.append(edge)
			else:
				right_node_edges.append(edge)

		edges_to_remove.extend(border_edges2)
		edges_to_remove.extend(right_chunk_edges)

		'''
		right_chunk_edges = set(searched_nodes[-2].get_edges())
		right_chunk_edges.difference_update(right_edges)
		'''
		searched_nodes[-2]._edges = right_chunk_edges
		searched_nodes[-2]._edges_sorted = False

		debug_temp2 = searched_nodes[-2].asm.high()

		right_trim_dist = searched_nodes[-2].asm.high() - region_hi 
		right_ref_CL = searched_nodes[-2].ref.trim_hi(right_trim_dist)
		right_asm_CL = searched_nodes[-2].asm.trim_hi(right_trim_dist)
		assert( right_asm_CL.low() == region_hi )
		assert( right_asm_CL.high() == debug_temp2 )
		debug_temp = len(right_asm_CL)
		right_asm_CL.left = region_lo
		right_asm_CL.right = searched_nodes[-2].asm.high() - chunk_len
		assert( debug_temp == len(right_asm_CL) )
		right_asm_og_CL = searched_nodes[-2].asm_original.trim_hi(right_trim_dist)
		right_node = Node(-1, right_ref_CL, right_asm_CL, right_asm_og_CL, right_node_edges)
		for edge in right_node_edges:
			assert( edge.edge_low(searched_nodes[-2]) >= right_node.asm.low() )
			assert( edge.edge_high(searched_nodes[-2]) <= right_node.asm.high() )
		right_seq = searched_nodes[-2].seq[(right_trim_dist + 1):]
		right_node.seq = right_seq
		assert len(right_seq) == len(right_asm_CL)
		searched_nodes[-2].seq = searched_nodes[-2].seq[:(right_trim_dist + 1)]
		assert len(searched_nodes[1].seq) == len(searched_nodes[1].asm)


	edge_list.update(edges_to_remove)
	b_e_temp_set = set(bad_edges)
	b_e_temp_set.difference_update(edge_list) #anything in edge_list that's also in b_e will be removed from b_e
	bad_edges = list(b_e_temp_set)

	bad_edges.sort( key = lambda e : e.edge_low(e.node1) ) #adding to make this stable/make the overall algorithm non-deterministic:
	#I have a bug that shows up many iterations into the program, at seemingly random times, so it's difficult to debug manually with prints
	#or with pdb. adding this to try to make the problem completely reproducible so I can hopefully finally track this down


	############### setting/updating prev and next for all affected nodes ################

	#logic explained in pic from 7/5/18

	if left_node_exists:
		leftmost = left_node
	elif right_node_exists:
		leftmost = right_node
	else:
		leftmost = searched_nodes[-1]

	if right_node_exists:
		rightmost = right_node
	elif left_node_exists:
		rightmost = left_node
	else:
		rightmost = searched_nodes[0]

	if leftmost is not None:
		leftmost.prev = searched_nodes[0]
	if rightmost is not None:
		rightmost.next = searched_nodes[-1]
	if searched_nodes[0] is not None:
		searched_nodes[0].next = leftmost
	if searched_nodes[-1] is not None:
		searched_nodes[-1].prev = rightmost
	if leftmost is not None and rightmost is not None and leftmost.asm.left < rightmost.asm.left:
		leftmost.next = rightmost
		rightmost.prev = leftmost
	searched_nodes[1].prev = insert_left
	if insert_left is not None:
		insert_left.next = searched_nodes[1]
	searched_nodes[-2].next = other_node
	other_node.prev = searched_nodes[-2]

	#update $scaffolds
	if searched_nodes[0] is None:
		if searched_nodes[-1] is None and not left_node_exists and not right_node_exists:
			#entire scaffold was placed within another
			del scaffolds[searched_nodes[1]]
		else:
			#remove old head/add new one
			for index, head in enumerate(scaffolds):
				if head is searched_nodes[1]:
					scaffolds[index] = leftmost
					break
			else:
				assert(1==2)
	if insert_left is None:
		#remove old head/add new one
		for index, head in enumerate(scaffolds):
			if head is other_node:
				scaffolds[index] = searched_nodes[1]
				break
		else:
			assert(1==2)



	#update edge endpoints (just nodes, not coordinates)
	if left_node_exists:
		left_node.new_edge_endpoints(searched_nodes[1])
	if right_node_exists:
		right_node.new_edge_endpoints(searched_nodes[-2])

	print("shift new chunk")
	chunk_node = searched_nodes[1]
	shift_dist = other_node.asm.low() - chunk_node.asm.low()
	new_name = other_node.asm.name

	debug_enum = 0
	while chunk_node is not other_node:
		#print(debug_enum)
		chunk_node.shift(shift_dist)
		chunk_node.asm.name = new_name
		for edge in chunk_node.get_edges():
			if edge.node1 is chunk_node:
				edge.asm1.name = new_name
			elif edge.node2 is chunk_node:
				edge.asm2.name = new_name
			else:
				assert 10 == 20
			assert edge.edge_low(chunk_node) >= chunk_node.asm.low()
			assert edge.edge_high(chunk_node) <= chunk_node.asm.high()
		chunk_node = chunk_node.next
		debug_enum += 1

	print("shift chunk's new right")
	shift_distance = (searched_nodes[-2].asm.high() - other_node.asm.low()) + 1
	new_right_node = other_node
	debug_enum2 = 0
	while new_right_node is not None:
		#print(debug_enum2)
		new_right_node.shift(shift_distance)
		for edge in new_right_node.get_edges():
			assert edge.edge_low(new_right_node) >= new_right_node.asm.low()
			assert edge.edge_high(new_right_node) <= new_right_node.asm.high()
		new_right_node = new_right_node.next
		debug_enum2 += 1

	if right_node_exists:
		if right_node.prev is not None:
			shift_dist = (right_node.prev.asm.high() + 1) - right_node.asm.low()
		else:
			shift_dist = -(right_node.asm.low() - 1) #want right node to start at 1 since its a scaffold head in this case
		right_node.shift_edges(shift_dist)
		for edge in right_node.get_edges():
			assert edge.edge_low(right_node) >= right_node.asm.low()
			assert edge.edge_high(right_node) <= right_node.asm.high()
	
	if rightmost is not None:
		pre_iter = rightmost.next
	elif searched_nodes[-1] is not None:
		pre_iter = searched_nodes[-1]
	else:
		pre_iter = None

	if pre_iter is not None:
		if pre_iter.prev is not None:
			shift_dist = (pre_iter.prev.asm.high() + 1) - pre_iter.asm.low()
		else:
			shift_dist = -(pre_iter.asm.low() - 1) #want right node to start at 1 since its a scaffold head in this case

	print("shift chunk's old right")
	while pre_iter is not None:
		pre_iter.shift(shift_dist)
		for edge in pre_iter.get_edges():
			assert edge.edge_low(pre_iter) >= pre_iter.asm.low()
			assert edge.edge_high(pre_iter) <= pre_iter.asm.high()
		pre_iter = pre_iter.next


	for head in scaffolds:
		while head is not None:
			for edge in head.get_edges():
				assert edge.edge_low(head) >= head.asm.low()
				assert edge.edge_high(head) <= head.asm.high()
			head = head.next







