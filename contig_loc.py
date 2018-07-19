class ContigLocation:

	def __init__(self, name, left, right):

		self.name = name
		self.left = left
		self.right = right

	def low(self):
		return min(self.left, self.right)

	def high(self):
		return max(self.left, self.right)

	def __str__(self):
		return self.name + "\t" + str(self.left) + "\t" + str(self.right)

	def shift(self, num):

		'''
		if self.rev():
			self.left = self.left - num
			self.right = self.right - num
		else:
			self.left = self.left + num
			self.right = self.right + num
		'''
		self.left = self.left + num
		self.right = self.right + num

	def midpoint(self):
		return (self.left + self.right) / 2

	def same_contig(self, otherCL):
		return self.name == otherCL.name

	def get_coords(self):
		return (self.left, self.right)

	def rev(self):
		return self.left > self.right

	#TODO newest algorithm shouldn't need this method; review and remove after completion if this is the case
	def trim(self, left_num, right_num):
		assert(6==7)
		temp = self.trim_left(left_num)
		return temp.trim_right(right_num)

	def trim_lo(self, num):

		#since both CLs shouldn't have the same start coord, and num refers to a point inclusively in a chunk region,
		#decrementing num prior to the construction of trimmed_from_left allows num to be used as an inclusive stop;
		#it is then immediately incremented back to its initial value so it can be used as an inclusive start
		if self.rev():
			num -= 1
			trimmed_from_left = ContigLocation(self.name, (self.right + num), self.right )
			num += 1
			self.right = self.right + num
			return trimmed_from_left

		else:
			num -= 1
			trimmed_from_left = ContigLocation(self.name, self.left, (self.left + num) )
			num += 1
			self.left = self.left + num
			return trimmed_from_left

	def trim_hi(self, num):

		#since both CLs shouldn't have the same start coord, and num refers to a point inclusively in a chunk region,
		#decrementing num prior to the construction of trimmed_from_right allows num to be used as an inclusive start;
		#it is then immediately incremented back to its initial value so it can be used as an inclusive stop
		if self.rev():
			num -= 1
			trimmed_from_right = ContigLocation(self.name, self.left, (self.left - num) )
			num += 1
			self.left = self.left - num
			return trimmed_from_right

		else:
			num -= 1
			trimmed_from_right = ContigLocation(self.name, (self.right - num), self.right )
			num += 1
			self.right = self.right - num
			return trimmed_from_right



	#TODO deprecated (maybe): the following 3 trim_no_change methods are for use with the oldest, most conservative version
	#of the chunk region algorithm, implemented in main.py at the time of writing. once the new algorithm is finished, the current
	#main.py and the following method definitions may be discarded
	def trim_no_change(self, left_num, right_num):
		temp = self.trim_left(left_num)
		return temp.trim_right(right_num)

	#TODO more accurately, these should be called trim_low and trim_high... aren't bidirectional alignments fun
	def trim_left_no_change(self, num):
		if self.rev():
			return ContigLocation(self.name, self.left, self.right + num)
		else:
			return ContigLocation(self.name, self.left + num, self.right)

	def trim_right_no_change(self, num):
		if self.rev():
			return ContigLocation(self.name, self.left - num, self.right)
		else:
			return ContigLocation(self.name, self.left, self.right - num)

	#TODO not sure if the +1 should be here- this means some portions of main.py can use the result directly, but some must subtract first
	#find all uses and determine the more common usage
	def __len__(self):
		return abs(self.right - self.left) + 1




