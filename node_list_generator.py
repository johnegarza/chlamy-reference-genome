from collections import defaultdict
import sys, os, csv, re, fnmatch
from block_node import Node
from fosmid_edge import Edge
import pickle

if len(sys.argv) < 3:
	sys.exit("Usage: %s block_list_tab_delimited indexed_fosmid_pairs" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])
if not os.path.exists(sys.argv[2]):
	sys.exit("Error: File '%s' not found" % sys.argv[2])


alignment_file = sys.argv[1] #tab_delim_results.tsv
fosmid_pairs = sys.argv[2]

contigs = [] #contains head node for each contig
line_indexed_nodes = [] #to retrieve node at line n, call line_indexed_nodes[n-1]
edges = [] #adding to help with active development

with open(alignment_file) as a_f:
	alignment_data = csv.reader(a_f, delimiter="\t")

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

for edge in edges:
	if ( edge.weight == -10 ):
		print("Edge between nodes in different scaffolds found")
		edge.this_node.print_surround_nodes()
		edge.other_node.print_surround_nodes()

#for some_node in line_indexed_nodes:
#	some_node.printn()
#	for some_edge in some_node.edges:
#		some_edge.printe(some_node)
#	print("")

#with open('linked_nodes.pkl', 'wb') as f1:
#	pickle.dump(contigs, f1, pickle.HIGHEST_PROTOCOL)
#with open('array_nodes.pkl', 'wb') as f2:
#	pickle.dump(line_indexed_nodes, f, pickle.HIGHEST_PROTOCOL)





