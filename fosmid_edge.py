class Edge:

	def __init__(self, this_node, other_node, t_r_s, t_r_e, t_a_s, t_a_e, o_r_s, o_r_e, o_a_s, o_a_e):

		self.this_node = this_node
		self.other_node = other_node

		self.this_ref_start = t_r_s
		self.this_ref_end = t_r_e
		self.this_asm_start = t_a_s
		self.this_asm_end = t_a_e
		self.other_ref_start = o_r_s
		self.other_ref_end = o_r_e
		self.other_asm_start = o_a_s
		self.other_asm_end = o_a_e

		start_avg = (self.this_asm_start + self.this_asm_end)/2
		stop_avg = (self.other_asm_start + self.other_asm_end)/2
		self.length = abs( start_avg-stop_avg )

		if (self.this_node.asm_name != self.other_node.asm_name):
			self.weight = -10
		elif ( (abs(self.length) < 35000) or (abs(self.length) > 40000) ):
			self.weight = 5
		else:
			self.weight = 10


	#deprecated?
	def printe(self, curr_node):
		
		if(self.this_node.line_num == curr_node.line_num):
			print("match to " + str(self.other_node.line_num) + " with weight " + str(self.weight) + "(" + str(abs(self.length)) + ")")
		else:
			print("match to " + str(self.this_node.line_num) + " with weight " + str(self.weight) + "(" + str(abs(self.length)) + ")")

