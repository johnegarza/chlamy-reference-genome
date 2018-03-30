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

#	def __str__(self):

	def tests(self):
		for edge in self.edges:

			if (self == edge.this_node):
				assert(self.ref_start <= edge.this_ref_start)
				assert(self.ref_stop >= edge.this_ref_end)
				#assembly alignments may be in reverse order
				block_asm_start = min(self.asm_start, self.asm_stop)
				block_asm_stop = max(self.asm_start, self.asm_stop)
				edge_asm_start = min(edge.this_asm_start, edge.this_asm_end)
				edge_asm_stop = max(edge.this_asm_start, edge.this_asm_end)
				assert(block_asm_start <= edge_asm_start)
				assert(block_asm_stop >= edge_asm_stop)
			else:
				assert(self.ref_start <= edge.other_ref_start)
				assert(self.ref_stop >= edge.other_ref_end)
				#assembly alignments may be in reverse order
				block_asm_start = min(self.asm_start, self.asm_stop)
				block_asm_stop = max(self.asm_start, self.asm_stop)
				edge_asm_start = min(edge.other_asm_start, edge.other_asm_end)
				edge_asm_stop = max(edge.other_asm_start, edge.other_asm_end)
				assert(block_asm_start <= edge_asm_start)
				assert(block_asm_stop >= edge_asm_stop)
				


