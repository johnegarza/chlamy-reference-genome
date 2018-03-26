#this file defines an object to support the one to one alignment data structure being built in delta_parser.py

class RefMapper:

	def __init__(self):

		#each element in the following 2 lists is a list of tuples; the overall list represents a block
		#of the reference and query that nucmer aligned, with the tuples specifying the exact 1-1 alignment
		self.ref_coords = []
		self.query_coords = []
		#this holds the name of the query scaffold/chromosome corresponding to a particular block
		self.query_names = []

		#TODO find a more effecient solution, or add in an option to disable this at the cost of increased computation time
		#reduces sort() time at the cost of memory; may be rebalanced in the future
		self.sort_amortizer = []

	def update(self, ref_coord_list, query_coord_list, q_name):
		self.ref_coords.append(ref_coord_list)
		self.query_coords.append(query_coord_list)
		self.query_names.append(q_name)

		self.sort_amortizer.append( (ref_coord_list, query_coord_list, q_name) )

	def __str__(self):
		
		#self.sort()
		for index, name in enumerate(self.query_names):

			print("----------------------------" + str(index) + "----------------------")

			print(str(name))
			print(self.ref_coords[index])
			print(self.query_coords[index])
			return "" #here because __str__ is required to return a string

	#this method sorts the reference coordinates list in ascending order
	#query coordinates and names are reordered along with it, preserving the mappings
	def sort(self):
		
		#0th element of each tuple in the list, which is a list of ref coords
		#0th element of a ref coords list is a tuple indicating the first alignment block
		#0th element of that tuple is the start coordinate of this entire list; since these are always
		#increaing, sorting by this value allows us to reconstruct the other 3 instance variables in the
		#sorted order described by the comments above this method's declaration
		self.sort_amortizer.sort(key = lambda triple : triple[0][0][0])

		#reset these, then rebuild in the for loop		
		self.ref_coords = []
		self.query_coords = []
		self.query_names = []
		for index, triple in enumerate(self.sort_amortizer):
			self.ref_coords.append(triple[0])
			self.query_coords.append(triple[1])
			self.query_names.append(triple[2])




