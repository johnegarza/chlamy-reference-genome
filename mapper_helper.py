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
		#actually, may be able to refactor the entire class to only need this list, and none of the 3 preceding ones
		self.sort_amortizer = []

	def update(self, ref_coord_list, query_coord_list, q_name):
		self.ref_coords.append(ref_coord_list)
		self.query_coords.append(query_coord_list)
		self.query_names.append(q_name)

		self.sort_amortizer.append( (ref_coord_list, query_coord_list, q_name) )

	def __str__(self):
		
		#self.sort()
		for index, name in enumerate(self.query_names):

			print(str(name))
			print(self.ref_coords[index])
			print(self.query_coords[index])
		return "" #here because __str__ is required to return a string

	def print_aligns(self, ref_name, line_num):

		for q_index, q_name in enumerate(self.query_names):

			#a for alignment?
			#can enumerate through ref or query coords, just chose ref
			for a_index, a_ref_block in enumerate(self.ref_coords[q_index]):

				join_list = []

				a_query_block = self.query_coords[q_index][a_index]

				join_list.append(ref_name)
				join_list.append(str(a_ref_block[0]))
				join_list.append(str(a_ref_block[1]))
				join_list.append(q_name)
				join_list.append(str(a_query_block[0]))
				join_list.append(str(a_query_block[1]))
				join_list.append(str(line_num))

				print("\t".join(join_list))

				line_num += 1

		return line_num


			#print(str(name) + "\t" + str(self.ref_coords[index]) + "\t" + str(self.query_coords[index]))

	#this method sorts the reference coordinates list in ascending order
	#query coordinates and names are reordered along with it, preserving the mappings
	def sort(self):
		
		#0th element of each tuple in the list is a list of ref coords
		#0th element of a ref coords list is a tuple indicating the first alignment block
		#0th element of that tuple is the start coordinate of this entire list; since these are always
		#increasing, sorting by this value allows us to reconstruct the other 3 instance variables in the
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

	#returns a 5-tuple: boolean for whether or not the input coordinates were successfully mapped, 
	#input left coord, input right coord, mapped scaffold name, mapped left coord, mapped right coord
	def map(self, left_coord, right_coord, line_num, num):

		list_index = -1
		for block_index, block in enumerate(self.ref_coords):

			#TODO remove int cast in final version of script, after all testing is complete
			#casting in order to make this explicitly fail while still under development, since in python 2.x,
			#comparing a str and int does not throw an error- it evaluates with no error and moves on
			if int(block[0][0]) <= int(left_coord):
				#list_index = block_index
				last_index = block[len(block)-1][1]
				if int(last_index) >= int(right_coord):
					#the fosmid end aligns to a continuous alignment block
					list_index = block_index
					break

		#TODO next 3 return statements are essentially filters; need to decide if these should silently
		#fail or not- data on which ends could not be used may be useful
		if (list_index == -1):
			#it was never updated
#			print str(line_num) + "." + str(num) + " has no match"
			return (False, -1, -1, "", -1, -1)

		
		left_tuple_index = -1
		right_tuple_index = -1
		for tuple_index, coord_tuple in enumerate(self.ref_coords[list_index]):
			if ((coord_tuple[0] <= left_coord) and (coord_tuple[1] >= left_coord)):
				left_tuple_index = tuple_index
			if ((coord_tuple[0] <= right_coord) and (coord_tuple[1] >= right_coord)):
				right_tuple_index = tuple_index

		if ( (left_tuple_index == -1) or (right_tuple_index == -1) or (left_tuple_index > right_tuple_index) ):
			
			#print( "left coord: " + str(left_coord) + " | right coord: " + str(right_coord) )
			#print( "left index: " + str(left_tuple_index) + " | right index: " + str(right_tuple_index) )
			#print(self.ref_coords[list_index])
			#print(self.query_coords[list_index])
			return (False, -1, -1, "", -1, -1)

		left_dist = left_coord - int(self.ref_coords[list_index][left_tuple_index][0])
		right_dist = right_coord - int(self.ref_coords[list_index][right_tuple_index][0])

		#check if the query is aligned forward or reverse
		if ( int(self.query_coords[list_index][left_tuple_index][0]) < int(self.query_coords[list_index][left_tuple_index][0]) ):
			newLeft = int(self.query_coords[list_index][left_tuple_index][0]) + left_dist
			newRight = int(self.query_coords[list_index][right_tuple_index][0]) + right_dist
		else:
			newLeft = int(self.query_coords[list_index][left_tuple_index][0]) - left_dist
			newRight = int(self.query_coords[list_index][right_tuple_index][0]) - right_dist


		old_dist = right_coord - left_coord
		new_dist = abs(newLeft - newRight)

		if ( abs(new_dist - old_dist) > 20 ):
			return (False, -1, -1, "", -1, -1)
			#print( str(left_coord) + "\t" + str(right_coord) + "\t" + str(old_dist) )
			#print( str(newLeft) + "\t" + str(newRight) + "\t" + str(new_dist) )

		return (True, left_coord, right_coord, self.query_names[list_index], newLeft, newRight)



