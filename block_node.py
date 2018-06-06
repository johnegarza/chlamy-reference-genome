from contig_loc import ContigLocation

class Node:

	def __init__(self, line_num, ref_name, ref_start, ref_stop, asm_name, asm_start, asm_stop):

		#line number in the nucmer alignment output- to keep track of where this came from and for debugging later
		self.line_num = int(line_num) 
		
		self.ref = ContigLocation(ref_name, ref_start, ref_stop)
		self.asm_original = ContigLocation(asm_name, asm_start, asm_stop)
		self.asm = ContigLocation(asm_name, asm_start, asm_stop)

		#the edges between this node and any others
		self._edges = []
		self._edges_sorted = True #used to help amortize sorting


		#since the final product should be a collection of nodes with at most one successor and one
		#predecessor, it works well as a linked list

		#TODO add these as parameters in the arguments list, with default values None
		self.prev = None
		self.next = None

	def add_edge(self, edge):
		self._edges.append(edge)
		self._edges_sorted = False

	def remove_edge(self, edge):
		self._edges.remove(edge)

	def get_edges(self):
		return self._edges

	#TODO double check to make sure this is sorting the way I expect it to
	#TODO quick familiar implementation below; this can be optimized further according to
	#https://stackoverflow.com/questions/403421/how-to-sort-a-list-of-objects-based-on-an-attribute-of-the-objects
	def get_sorted_edges(self):
		if (not _edges_sorted):
			self._edges.sort( key = lambda e : e.edge_start(self) )
			_edges_sorted = True
		return self._edges

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

		for edge in self._edges:

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
