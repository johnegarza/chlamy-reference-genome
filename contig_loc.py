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

		if self.rev():
			self.left = self.left - num
			self.right = self.right - num
		else:
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

	def trim(self, left_num, right_num):
		temp = self.trim_left(left_num)
		return temp.trim_right(right_num)

	def trim_left(self, num):
		if self.rev():
			return ContigLocation(self.name, self.left - num, self.right)
		else:
			return ContigLocation(self.name, self.left + num, self.right)

	def trim_right(self, num):
		if self.rev():
			return ContigLocation(self.name, self.left, self.right + num)
		else:
			return ContigLocation(self.name, self.left, self.right - num)

	def __len__(self):
		return abs(self.right - self.left)




