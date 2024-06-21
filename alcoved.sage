R = RootSystem(['A',2])

def AlcovedPolytope:
	def __init__(self, R, BoundaryParameters):
		self.root_system = R
		self.boundaries = BoundaryParameters
		self.fundamental_coweights = R.coweight_lattice().basis()
		self.rho = sum(self.fundamental_coweights)
		self.dual_Coxeter_no = R.root_lattice().highest_root().dot_product(self.rho)+1
		
	def is_member(self, coweight):
		for key, value in self.boundaries:
			foo = coweight.dot_product(key)
			if foo < value[0] or foo > value[1]:
				return False
		return True

	def cover(self, coweight):
		if !self.is_member(coweight):
			raise Exception('coweight is not in the alcoved polytope')
		foo = self.fundamental_coweights.list().append(self.rho)
		return sum([i for i in foo if self.is_member(coweight-i)])