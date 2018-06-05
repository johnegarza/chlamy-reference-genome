class ContigLocation:

	def __init(self, name, left, right):

		self.name = name
		self.left = left
		self.right = right

		self.low = min(left, right)
		self.high = max(left, right)
