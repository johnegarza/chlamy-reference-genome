from contig_loc import ContigLocation

class Edge:

	def __init__(self, this_node, other_node, t_r_s, t_r_e, t_a_s, t_a_e, o_r_s, o_r_e, o_a_s, o_a_e):

		#Endpoint 1
		self.ref1 = ContigLocation("xxx", t_r_s, t_r_e)
		self.asm1_original = ContigLocation("xxx", t_a_s, t_a_e)
		self.asm1 = ContigLocation("xxx", t_a_s, t_a_e)
		self.node1 = this_node

		#Endpoint 2
		self.ref2 = ContigLocation("xxx", o_r_s, o_r_e)
		self.asm2_original = ContigLocation("xxx", o_a_s, o_a_e)
		self.asm2 = ContigLocation("xxx", o_a_s, o_a_e)
		self.node2 = other_node

		self.length = abs ( self.asm1.midpoint() - self.asm2.midpoint() )

		if ( not self.asm1.same_contig(self.asm2) ):
			self.weight = -10
		elif ( (self.length < 35000) or (self.length > 40000) ):
			self.weight = 5
		else:
			self.weight = 10
		'''
		#temp testing
		assert t_a_s <= t_a_e, str(t_a_s) + " " + str(t_a_e)
		assert o_a_s <= o_a_e

	def edge_start(self, node):
		if node is self.node1:
		elif node is self.node2:
		else:
			assert (1==2)
		'''

	def edge_low(self, node):
		if node is self.node1:
			return self.asm1.low()
		elif node is self.node2:
			return self.asm2.low()
		else:
			assert(1==2)

	def edge_high(self, node):
		if node is self.node1:
			return self.asm1.high()
		elif node is self.node2:
			return self.asm2.high()
		else:
			assert(1==2)

	#edge_start/stop may be deprecated since the introduction of edge_low/high
	def edge_start(self, node):
		if node is self.node1:
			return self.asm1.left
		elif node is self.node2:
			return self.asm2.left
		else:
			assert(1==2)

	def edge_stop(self, node):
		if node is self.node1:
			return self.asm1.right
		elif node is self.node2:
			return self.asm2.right
		else:
			assert(1==2)

	def opposite_node(self, node):
		if node is self.node1:
			return self.node2
		elif node is self.node2:
			return self.node1
		else:
			assert(1==2)

	def other_node_info(self, node):
		if(node is self.node1):
			return str(self.node2)
		elif(node is self.node2):
			return str(self.node1)
		else: #shouldn't happen
			assert(1==2)

	#called from tests() function in block_node.py
	def node_info(self, node):
		ans = []
		if(node is self.this_node):
			ans.append(str(self.this_ref_start))
			ans.append(str(self.this_ref_end))
			ans.append(str(self.this_asm_start))
			ans.append(str(self.this_asm_end))
			return "\t".join(ans)
		elif(node is self.other_node):
			ans.append(str(self.other_ref_start))
			ans.append(str(self.other_ref_end))
			ans.append(str(self.other_asm_start))
			ans.append(str(self.other_asm_end))
			return "\t".join(ans)

		else:
			#this shouldn't happen, but for safety...
			return "Bad node passed to edge"

	def __str__(self):
		return str(self.node1.asm) + "\t" + (self.node2.asm)


