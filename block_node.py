from contig_loc import ContigLocation

class Node:

	def __init__(self, line_num, ref_name, ref_start, ref_stop, asm_name, asm_start, asm_stop):

		#line number in the nucmer alignment output- to keep track of where this came from and for debugging later
		self.line_num = int(line_num) 
		
		self.ref = ContigLocation(ref_name, ref_start, ref_stop)
		self.asm_original = ContigLocation(asm_name, asm_start, asm_stop)
		self.asm = ContigLocation(asm_name, asm_start, asm_stop)

		#the edges between this node and any others
		self.edges = []

		#since the final product should be a collection of nodes with at most one successor and one
		#predecessor, it works well as a linked list
		self.prev = None
		self.next = None


	def printn(self):
		print(str(self.line_num))

	def __str__(self):
		ans = []
		ans.append(str(self.line_num))
		ans.append(str(self.ref))
		ans.append(str(self.asm))
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

			#assembly alignments may be in reverse order
			block_asm_start = self.asm.low()
			block_asm_stop = self.asm.high()

			if (self is edge.node1):

				edge_asm_start = edge.asm1.low()
				edge_asm_stop = edge.asm1.high()

				if not(block_asm_start <= edge_asm_start):
					debug_print = True
				if not(block_asm_stop >= edge_asm_stop):
					debug_print = True

			elif(self is edge.node2):

				edge_asm_start = edge.asm2.low()
				edge_asm_stop = edge.asm2.high()

				if not(block_asm_start <= edge_asm_start):
					debug_print = True
				if not(block_asm_stop >= edge_asm_stop):
					debug_print = True
			else:
				assert(1==2) #make this fail noticeably- this case shouldn't ever happen during normal execution

			if(debug_print):
				print(edge.node_info(self))
