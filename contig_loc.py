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
