class ContigLocation:

	def __init__(self, name, left, right):

		self.name = name
		self.left = left
		self.right = right

		self.low = min(left, right)
		self.high = max(left, right)

	def shift(self, num):

		self.left = self.left + num
		self.right = self.right + num
		self.low = self.low + num
		self.high = self.high + num
