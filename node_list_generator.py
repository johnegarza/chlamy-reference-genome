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
#	python node_list_generator.py htcf_data/node_list.tsv delta_mapped_ends.txt

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
edges = [] #adding to help with active development

with open(alignment_file) as a_f:
	alignment_data = csv.reader(a_f, delimiter="\t")

	#TODO refactor to use Contig (upcoming) data structure, so don't need sentinel nodes

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
		node1.add_edge(edge)
		node2.add_edge(edge)
		edges.append(edge) #TODO if no proper use for this, remove; will just lead to memory leaks, as this keeps edges deleted later on still alive due to the reference

timer = len(line_indexed_nodes)
for num, node in enumerate(line_indexed_nodes):

	print( str(num+1) + "/" + str(timer) )

	bad_edges = set()
	chunk_seeds = []

	#find all other assembly scaffolds this node may have edges to; will choose one edges for each different scaffold, then use that as a "seed" to build up 
	#a "chunk", a portion of the node that likely belongs on that other scaffold and not the current one
	for edge in node.get_edges():
		try:
			a_name = edge.opposite_node(node).asm.name
		except AttributeError:
			print(edge)

		if ( (edge.weight == -10) and (a_name not in bad_edges) ):
			bad_edges.add(a_name)
			chunk_seeds.append(edge)

	###begin new dev

	#TODO how well does this work for ref-asm reversed alignments

	this_asm_name = node.asm.name
#	test_num = 0
	for seed in chunk_seeds:
#		test_num += 1

		#this is necessary because after the first iteration of the loop, the original "node" now exists; but we still have the edges, which have 
		#been updated and can be used to recover the node containing the end we care about
		node = seed.node_by_name(this_asm_name)

		#STEP 1: grab all edges that belong in a chunk and define chunk boundaries

		#without slice, this is just a reference to the node edge list, so when things are removed from there they're removed from ordered edges too
		#TODO make sure this problem isn't happening anywhere else- see other todo labeled "search me"
		ordered_edges = node.get_sorted_edges()[:] 
		other_node = seed.opposite_node(node)

		chunk_edges = set()
		index = -1
		for loc, search in enumerate(ordered_edges):
			if search is seed:
				index = loc
				break
		chunk_start = seed.edge_low(node)
		chunk_stop = seed.edge_high(node)

		#TODO this search can be optimized- currently building for speed in order to demo
		#TODO is len inside a loop condition recalculated every time? if so could be made more efficient, but is it possible for the list

		left_search = index - 1
		while(left_search >= 0):
			curr_edge = ordered_edges[left_search]
			if curr_edge.opposite_node(node) is other_node:
				chunk_edges.add(curr_edge)
			else:
				break #chunk is a contiguous section of a node containing edges that all end in the same other node, which is on a different contig
				#thus, encountering an edge that does not end in that particular different contig indicates a chunk boundary
			chunk_start = min( chunk_start, curr_edge.edge_low(node) )
			chunk_stop = max( chunk_stop, curr_edge.edge_high(node) )
			left_search -= 1

		right_search = index + 1
		#to be modified (anticipating future parallelization of portions of the code)
		length = len(ordered_edges)
		while ( right_search < length):
			curr_edge = ordered_edges[right_search]
			if curr_edge.opposite_node(node) is other_node:
				chunk_edges.add(curr_edge)
			else:
				break #chunk is a contiguous section of a node containing edges that all end in the same other node, which is on a different contig
				#thus, encountering an edge that does not end in that particular different contig indicates a chunk boundary
			chunk_start = min( chunk_start, curr_edge.edge_low(node) )
			chunk_stop = max( chunk_stop, curr_edge.edge_high(node) )
			right_search += 1

		chunk_edges.add(seed)

		#right_search and left_search are exclusive coordinates relative to the edges placed in chunk_edges from ordered_edges
		assert( len(chunk_edges) == (right_search - left_search - 1) )

		#STEP 2 separate the remaining edges to make sure none accidentally have an endpoint cut in half when the chunk is removed
		#for reference, see pic from 6/6/18 ~12:15 PM TODO download and include in repo
		
		right_edges = [] #this ensures proper behavior in the event that right_search = ordered_edges and the loop doesn't run
		while ( right_search < len(ordered_edges) ):
			curr_edge = ordered_edges[right_search]
			if curr_edge.edge_low(node) < chunk_stop :
				node.remove_edge(curr_edge)
				curr_edge.opposite_node(node).remove_edge(curr_edge)
#				print("remove edge " + str(curr_edge) )
				#TODO
				#is this sufficient? should curr_edge be explicitly deleted as well?
				#should the removal of an edge be noted? this may be dangerous
			else:
				right_edges = ordered_edges[right_search:]
				break
			right_search += 1

		left_edges = ordered_edges[:left_search] #TODO changed from left_edges
		while ( left_search >= 0 ):

			

			curr_edge = ordered_edges[left_search]
			if curr_edge.edge_high(node) > chunk_start :
				#TODO "search me" before slice, the following remove was removing from ordered_edges too
				node.remove_edge(curr_edge)
				curr_edge.opposite_node(node).remove_edge(curr_edge)
#				print("remove edge " + str(curr_edge) )
				#TODO
				#is this sufficient? should curr_edge be explicitly deleted as well?
				#should the removal of an edge be noted? this may be dangerous
				
				#TODO should I be doing this
				try:
					left_edges.remove(curr_edge)
				except ValueError:
					pass

			#no else- in this while case, we must check until the end, due to the ordering being based on the
			#left/start coord (see pic for reference- must remove edges 3 and 4, and can't know how many
			#have stop coords after chunk start from the start coord alone)

			left_search -= 1

		#safety check while developing TODO investigate why this isn't working
		#assert len(node.get_edges()) == ( len(chunk_edges) + len(right_edges) + len(left_edges) )
		if len(node.get_edges()) != ( len(chunk_edges) + len(right_edges) + len(left_edges) ):
			print(len(node.get_edges()))
			print(len(chunk_edges))
			print(len(right_edges))
			print(len(left_edges))

			moved_edges = []
			moved_edges.extend(chunk_edges)
			moved_edges.extend(right_edges)
			moved_edges.extend(left_edges)
			lost_edges = list(set(node.get_edges()).difference(moved_edges))
			for thing in lost_edges:
				print(thing)
		
		#STEP 3 construct new nodes

		#TODO this approach simply inserts a pulled-out bad node to the left of the node that it matches- real algorithm will be much more complex
		
		node_len = (chunk_stop - chunk_start) + 1
		#another = other.prev

#		print(test_num)
#		print(node.asm)
#		print(type(node.asm))
		left_dist = chunk_start - node.asm.left #using this gives inclusive coords
		right_dist = node.asm.right - chunk_stop

		#TODO refactor with new ContigLocation trim functions
		#TODO ref CL changes only work if nodes are a one-to-one mapping; is this accurate?
		if node.ref.rev():
			chunk_ref_CL = ContigLocation(node.ref.name, node.ref.left - left_dist, node.ref.right + right_dist)
		else:
			chunk_ref_CL = ContigLocation(node.ref.name, node.ref.left + left_dist, node.ref.right - right_dist)

		chunk_asm_CL = ContigLocation(other_node.asm.name, other_node.asm.left, other_node.asm.right + (node_len - 1) )

		chunk_node = Node(-1, chunk_ref_CL, chunk_asm_CL, node.asm_original, chunk_edges)

		left_ref_CL = node.ref.trim_left(left_dist - 1) #-1 because otherwise this and prev node would start at the exact same coord; this CL should have exclusive coords
		left_asm_CL = ContigLocation(node.asm.name, node.asm.left, chunk_start - 1)
		left_node = Node(-1, left_ref_CL, left_asm_CL, node.asm_original, left_edges)

		right_ref_CL = node.ref.trim_right(right_dist - 1)
		right_asm_CL = ContigLocation(node.asm.name, chunk_start, node.asm.right - node_len)
		right_node = Node(-1, right_ref_CL, right_asm_CL, node.asm_original, right_edges)

		#STEP 4 insert new nodes, including updating references for the surrounding nodes

		left_node.next = right_node
		right_node.prev = left_node

		#node could be the head of a contig
		if node.prev is not None:
			node.prev.next = left_node
			left_node.prev = node.prev
		else:
			for index, head in enumerate(contigs):
				if head is node:
					contigs[index] = left_node


		#node could be the tail of a contig
		if node.next is not None:
			node.next.prev = right_node
			right_node.next = node.next



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

		#STEP 5 update edge endpoints to point to their new nodes and clear original node so that no refs are left and garbage collector will free its memory
		left_node.new_edge_endpoints(node)
		right_node.new_edge_endpoints(node)
		chunk_node.new_edge_endpoints(node)
		node.clear()
#		del node #TODO evaluate if this is necessary when the design is more complete


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

		break


	###end new dev

	'''
	for seed in chunk_seeds:

		other = seed.opposite_node(node)
		a_name = other.asm.name

		new_start, new_stop = other.asm.get_coords()
		left_bound = new_stop
		right_bound = new_start
		new_edges = set() #to avoid accidental duplications, this is not a list

		#build up a region (chunk) that groups all of the edges a to particular, different contig
		#this will be pulled out, formed into a new node, and placed elsewhere
		for edge in node.get_edges():
			temp_start, temp_stop = edge.opposite_node(node).asm.get_coords()
			if( (edge.weight == -10) and (a_name == edge.opposite_node(node).asm.name ) and (temp_start < left_bound) ):
				new_start = temp_start
				new_edges.add(edge)
			if( (edge.weight == -10) and (a_name == edge.opposite_node(node).asm.name ) and (temp_stop > right_bound) ):
				new_stop = temp_stop
				new_edges.add(edge)
		new_edges.add(seed)

		#create the new node
		move_record = seed.node1.asm.name + ">" + seed.node2.asm.name
		#move_record = node.asm_name + ">" + other.asm_name #record where it was and where it was placed
		new_CL = ContigLocation(move_record, new_start, new_stop)
		new_node = Node(-1, node.ref, new_CL)

		#transfer ownership of edges
		for edge in new_edges:
			node.remove_edge(edge)
			new_node.add_edge(edge)
			if (edge.node1 is node):
				edge.node1 = new_node
			elif (edge.node2 is node):
				edge.node2 = new_node
			else:
				assert(3==4)

		#place node- naive approach for gradual buildup of algorithm
		#TODO TODO TODO TODO placement should be dynamic in later versions

		#had to move the following statement to its current location
		#since seed was being reassigned and thus failing the lookup (ret_other_node)
		#other = seed.ret_other_node(node)
		new_node.prev = other.prev
		new_node.next = other
		if(other.prev is not None): #TODO and is not sentinel node
			other.prev.next = new_node
		else:
			#creating a new contig head, so update the contigs list
			for index, contig in enumerate(contigs):
				if contig is other:
					contigs[index] = new_node #TODO forgot about sentinel nodes- this line does nothing

		other.prev = new_node
		
		chunk_length = abs(new_stop - new_start) + 1

		#adjust the stop coord of this node since we just pulled a chunk out
		node.asm_stop = node.asm_stop - chunk_length

		#adjust the node (and corresponding edge) assembly coords for all nodes in this contig following the altered one
		iterator = node.next
		while(iterator is not None):

			iterator.asm.shift(-chunk_length)

			for edge in iterator.get_edges():

				edge.opposite_node(iterator).asm.shift(-chunk_length)

			iterator = iterator.next

		#adjust the nodes that come after the new node
		iterator = new_node.next
		while(iterator is not None):
			iterator.asm.left = iterator.asm.left + chunk_length
			iterator.asm.right = iterator.asm.right + chunk_length

			for edge in iterator.get_edges():

				edge.opposite_node(iterator).asm.shift(chunk_length)

			iterator = iterator.next
		'''
'''		
for cnum, contig_head in enumerate(contigs):
	
	iterator = contig_head
	while(iterator is not None):
		print( str(cnum) + "\t" + str(iterator) )
		iterator = iterator.next
'''

for cnum, contig_head in enumerate(contigs):
	
	iterator = contig_head
	while(iterator is not None):
		for edge in iterator._edges:
			if edge.weight == -10:
				print("bad edge")
		iterator = iterator.next



