from collections import defaultdict
import sys, os, csv, re, fnmatch
from block_node import Node
from fosmid_edge import Edge
import pickle
import argparse
import copy

if len(sys.argv) < 3:
	sys.exit("Usage: %s block_list_tab_delimited indexed_fosmid_pairs" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])
if not os.path.exists(sys.argv[2]):
	sys.exit("Error: File '%s' not found" % sys.argv[2])


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
	dummy_node = Node(0, "dummmy", 0, 0, "dummy", 0, 0)
	curr_node = dummy_node

	for line_id, block in enumerate(alignment_data):

		#load in data from the current line
		ref_chr = block[0]
		ref_start = int(block[1])
		ref_stop = int(block[2])
		asm_scaf = block[3]
		asm_start = int(block[4])
		asm_stop = int(block[5])
		line_num = int(block[8])

		if(prev_scaf != asm_scaf): #end of an assembly contig

			#print("line: " + str(line_num) + " prev scaf: " + prev_scaf + " curr scaf: " + asm_scaf)

			#head/tail id both for use in debugging and to force creation of new node each time
			#without unique id, might just repeatedly use reference to same head/tail each time
			#so all contigs would share same head/tail and things would get weird (actually, only
			#the last contig would use them, and refs to all other contigs would be lost)

			tail_id = "tail after " + str(line_id)
			tail_node = Node(0, tail_id, 0, 0, tail_id, 0, 0)
			tail_node.prev = curr_node
			curr_node.next = tail_node

			head_id = "head before " + str(line_id + 1)
			head_node = Node(0, head_id, 0, 0, head_id, 0, 0)
			curr_node = Node(line_num, ref_chr, ref_start, ref_stop, asm_scaf, asm_start, asm_stop)
			head_node.next = curr_node
			curr_node.prev = head_node

			contigs.append(head_node)

			line_indexed_nodes.append(curr_node)

			prev_scaf = asm_scaf
			
			#TODO isn't this just an else case
			continue #work is done for this cycle, and thanks to python's scope (or lack thereof)
				 #we can just skip to the next iteration of the loop

		new_node = Node(line_num, ref_chr, ref_start, ref_stop, asm_scaf, asm_start, asm_stop)
		curr_node.next = new_node
		new_node.prev = curr_node
		curr_node = new_node

		line_indexed_nodes.append(curr_node)

		prev_scaf = asm_scaf

	#make the tail for the final contig
	tail_id = "tail after " + str(line_id)
	tail_node = Node(0, tail_id, 0, 0, tail_id, 0, 0)
	tail_node.prev = curr_node
	curr_node.next = tail_node

#for contig_head in contigs:
#	curr = contig_head
#	while(curr.next != None):
#		curr.printn()
#		curr = curr.next
#	curr.printn() # print the tails
#	print("switch")

#line_indexed_nodes[500].printn()
#line_indexed_nodes[1000].printn()
#line_indexed_nodes[1500].printn()
#line_indexed_nodes[2000].printn()

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

		#print(str(endpoint1) + "\t" + str(node1.line_num))

		edge = Edge(node1, node2, left_ref_start, left_ref_stop, left_asm_start, left_asm_stop, right_ref_start, right_ref_stop, right_asm_start, right_asm_stop)
		node1.edges.append(edge)
		node2.edges.append(edge)
		edges.append(edge)

#some basic tests
#for some_node in line_indexed_nodes:
#	some_node.tests()

#print(len(contigs))
#print(len(line_indexed_nodes))
#print(len(edges))
#for edge in edges:
#	if ( edge.weight == -10 ):
#		print("hit")
		#print("Edge between nodes in different scaffolds found")
		#edge.this_node.print_surround_nodes()
		#edge.other_node.print_surround_nodes()

#print(len(line_indexed_nodes))
#for node in line_indexed_nodes:
#	checker = True
#	for edge in node.edges:
#		if (edge.weight == -10):
#			checker = False
#	if(checker):
#		print("good node")


#for node in line_indexed_nodes:
#	checker = False
#	for edge in node.edges:
#		if (edge.weight == -10):
#			checker = True
#			break
#	if(checker):
#		print("bad node")

'''
for node in line_indexed_nodes:
	node_print = True
	for edge in node.edges:
		if (edge.weight == -10):
			if node_print:
				print(str(node))
				print("------")
				node_print = False
			print(edge.other_node_info(node))
	print("")
'''

#for some_node in line_indexed_nodes:
#	some_node.printn()
#	for some_edge in some_node.edges:
#		some_edge.printe(some_node)
#	print("")

#with open('linked_nodes.pkl', 'wb') as f1:
#	pickle.dump(contigs, f1, pickle.HIGHEST_PROTOCOL)
#with open('array_nodes.pkl', 'wb') as f2:
#	pickle.dump(line_indexed_nodes, f, pickle.HIGHEST_PROTOCOL)

sys.setrecursionlimit(10000) #TODO there has to be a better way to do this

alg_list = copy.deepcopy(line_indexed_nodes)

num_of_nodes = len(line_indexed_nodes)

for num, node in enumerate(line_indexed_nodes):

	#print( str(num+1) + " of " + str(num_of_nodes) )

	bad_edges = set()
	chunk_seeds = []

	for edge in node.edges:
		a_name = edge.other_node_asm_name(node)
		if ( (edge.weight == -10) and (a_name not in bad_edges) ):
#			print(edge)
			bad_edges.add(a_name)
			chunk_seeds.append(edge)

	for seed in chunk_seeds:
		a_name = seed.other_node_asm_name(node)
		other = seed.ret_other_node(node)

		new_start, new_stop = seed.other_node_asm_coords(node)
		left_bound = new_stop
		right_bound = new_start
		new_edges = set() #to avoid accidental duplications, this is not a list

		#build up a region (chunk) that groups all of the edges a to particular, different contig
		#this will be pulled out, formed into a new node, and placed elsewhere
		for edge in node.edges:
			temp_start, temp_stop = edge.other_node_asm_coords(node)
			if( (edge.weight == -10) and (a_name == edge.other_node_asm_name(node) ) and (temp_start < left_bound) ):
				new_start = temp_start
				new_edges.add(edge)
			if( (edge.weight == -10) and (a_name == edge.other_node_asm_name(node) ) and (temp_stop > right_bound) ):
				new_stop = temp_stop
				new_edges.add(edge)
		new_edges.add(seed)

		#create the new node
		move_record = seed.node1.asm.name + ">" + seed.node2.asm.name
		#move_record = node.asm_name + ">" + other.asm_name #record where it was and where it was placed
		new_node = Node(-1, node.ref.name, node.ref.left, node.ref.right, move_record, new_start, new_stop)

		#transfer ownership of edges
		for edge in new_edges:
			node.edges.remove(edge)
			new_node.edges.append(edge)
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
		
		chunk_length = (new_stop - new_start) + 1

		#adjust the stop coord of this node since we just pulled a chunk out
		node.asm_stop = node.asm_stop - chunk_length

		#adjust the node (and corresponding edge) assembly coords for all nodes in this contig following the altered one
		iterator = node.next
		while(iterator is not None):

#			iterator.asm_start = iterator.asm_start - chunk_length
#			iterator.asm_stop = iterator.asm_stop - chunk_length
			#TODO further modify to use the CL shift() method once updates are complete- same for next 2 for loops
			iterator.asm.left = iterator.asm.left - chunk_length
			iterator.asm.right = iterator.asm.right - chunk_length

			for edge in iterator.edges:

				if (edge.node1 is iterator):
#					edge.this_asm_start = edge.this_asm_start - chunk_length
#					edge.this_asm_end = edge.this_asm_end - chunk_length
					edge.asm1.left = edge.asm1.left - chunk_length
					edge.asm1.right = edge.asm1.right - chunk_length
				elif (edge.node2 is iterator):
					edge.asm2.left = edge.asm2.left - chunk_length
					edge.asm2.right = edge.asm2.right - chunk_length
				else:
					assert(1==2)

			iterator = iterator.next

		#adjust the nodes that come after the new node
		iterator = new_node.next
		while(iterator is not None):
			iterator.asm.left = iterator.asm.left + chunk_length
			iterator.asm.right = iterator.asm.right + chunk_length

			for edge in iterator.edges:

				if (edge.node1 is iterator):
					edge.asm1.left = edge.asm1.left + chunk_length
					edge.asm1.right = edge.asm1.right + chunk_length
				elif (edge.node2 is iterator):
					edge.asm2.left = edge.asm2.left + chunk_length
					edge.asm2.right = edge.asm2.right + chunk_length
				else:
					assert(1==2)

			iterator = iterator.next
		
for cnum, contig_head in enumerate(contigs):
	
	iterator = contig_head
	while(iterator is not None):
		print( str(cnum) + "\t" + str(iterator) )
		iterator = iterator.next




