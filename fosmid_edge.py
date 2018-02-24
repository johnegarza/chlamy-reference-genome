class Edge:

	def __init__(self, this_node, other_node):

		self.this_node = this_node
		self.other_node = other_node

		avg_endpoint1 = (this_node.asm_start + this_node.asm_stop)/2
		avg_endpoint2 = (other_node.asm_start + other_node.asm_stop)/2
		self.length = (avg_endpoint1 - avg_endpoint2)

		#debugging
		if(this_node.line_num == 673 or other_node.line_num == 673):
			print(str(avg_endpoint1) + "\t" + str(avg_endpoint2) + "\t" + str(self.length))


		if (self.this_node.asm_name != self.other_node.asm_name):
			self.weight = -10
		elif ( (abs(self.length) < 35000) or (abs(self.length) > 40000) ):
			self.weight = 5
		else:
			self.weight = 10


	def printe(self, curr_node):
		
		if(self.this_node.line_num == curr_node.line_num):
			print("match to " + str(self.other_node.line_num) + " with weight " + str(self.weight) + "(" + str(abs(self.length)) + ")")
		else:
			print("match to " + str(self.this_node.line_num) + " with weight " + str(self.weight) + "(" + str(abs(self.length)) + ")")

