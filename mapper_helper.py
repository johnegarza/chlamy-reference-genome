#this file defines an object to support the one to one alignment data structure being built in delta_parser.py

class RefMapper:

	def __init__(self):

		#each element in the following 2 lists is a list of tuples; the overall list represents a block
		#of the reference and query that nucmer aligned, with the tuples specifying the exact 1-1 alignment
		self.ref_coords = []
		self.query_coords = []
		#this holds the name of the query scaffold/chromosome corresponding to a particular block
		self.query_name = []

	def update(self, ref_coord_list, query_coord_list, q_name):
		self.ref_coords.append(ref_coord_list)
		self.query_coords.append(query_coord_list)
		self.query_name.append(q_name)

	def __str__(self):
		
		for index, name in enumerate(self.query_name):

			print(str(name))
			print(self.ref_coords[index])
			print(self.query_coords[index])
			return "" #here because __str__ is required to return a string
