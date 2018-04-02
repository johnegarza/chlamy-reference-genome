class Node:

	def __init__(self, line_num, ref_name, ref_start, ref_stop, asm_name, asm_start, asm_stop):

		#line number in the nucmer alignment output- to keep track of where this came from and for debugging later
		self.line_num = int(line_num) 
		self.ref_name = str(ref_name) #chromosome/scaffold name in phytozome reference
		self.ref_start = int(ref_start) #start and stop coordinates in reference
		self.ref_stop = int(ref_stop)
		self.asm_name = str(asm_name) #scaffold name in our assembly
		self.asm_start = int(asm_start) #start and stop coordinates in our assembly
		self.asm_stop = int(asm_stop)

		#the edges between this node and any others
		self.edges = [] #how? tuple? edge object?

		#since the final product should be a collection of nodes with at most one successor and one
		#predecessor, it works well as a linked list
		self.prev = None
		self.next = None


	def printn(self):
		print(str(self.line_num))

	def __str__(self):
		ans = []
		ans.append(str(self.line_num))
		ans.append(self.ref_name)
		ans.append(str(self.ref_start))
		ans.append(str(self.ref_stop))
		ans.append(self.asm_name)
		ans.append(str(self.asm_start))
		ans.append(str(self.asm_stop))
		return "\t".join(ans)

	def print_surround_nodes(self):
		prevN = "No prev node"
		nextN = "No next node"
		if self.prev is not None:
			prevN = str(self.prev)
		if self.next is not None:
			nextN = str(self.next)
		print(prevN)
		print(str(self))
		print(nextN)


	def tests(self):

		print(str(self))
		debug_print = False

		for edge in self.edges:

			debug_print = False

			if (self is edge.this_node):
#				if not(self.ref_start <= int(edge.this_ref_start)):
#					debug_print = True
#				if not(self.ref_stop >= int(edge.this_ref_end)):
#					debug_print = True
				#assembly alignments may be in reverse order
				block_asm_start = min(self.asm_start, self.asm_stop)
				block_asm_stop = max(self.asm_start, self.asm_stop)
				edge_asm_start = min(edge.this_asm_start, edge.this_asm_end)
				edge_asm_stop = max(edge.this_asm_start, edge.this_asm_end)
				if not(block_asm_start <= edge_asm_start):
					debug_print = True
				if not(block_asm_stop >= edge_asm_stop):
					debug_print = True

			elif(self is edge.other_node):
#				if not(self.ref_start <= int(edge.other_ref_start)):
#					debug_print = True
#				if not(self.ref_stop >= int(edge.other_ref_end)):
#					debug_print = True
				#assembly alignments may be in reverse order
				block_asm_start = min(self.asm_start, self.asm_stop)
				block_asm_stop = max(self.asm_start, self.asm_stop)
				edge_asm_start = min(edge.other_asm_start, edge.other_asm_end)
				edge_asm_stop = max(edge.other_asm_start, edge.other_asm_end)
				if not(block_asm_start <= edge_asm_start):
					debug_print = True
				if not(block_asm_stop >= edge_asm_stop):
					debug_print = True
			else:
				print("NOT POSSIBLE")

			if(debug_print):
				print(edge.node_info(self))
