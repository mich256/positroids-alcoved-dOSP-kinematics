def fundamental_coweights(R):
	return R.coweight_space().basis()

def positive_roots(R):
	theta = R.root_lattice().highest_root()
	coeff = theta.coefficients()
	n = len(R.index_set())
	sr = R.root_lattice().simple_roots()
	def helper(i):
		if i == n:
			return [j*sr[n] for j in (1..coeff[n-1])]
		foo = []
		for j in (1..coeff[i-1]):
			foo.append(j*sr[i])
		for k in helper(i+1):
			for j in (1..coeff[i-1]):
				foo.append(j*sr[i] + k)
		return foo
	bar = []
	for i in R.index_set():
		for j in helper(i):
			bar.append(j)
	return bar

class AlcovedPolytope:
	def __init__(self, R, BoundaryParameters):
		self.root_system = R
		self.boundaries = BoundaryParameters
		# the coweights are indexed from 1, not from 0
		self.quotient_basis = [fundamental_coweights(R)[i]/R.cartan_type().dual_coxeter_number() for i in R.index_set()]
		
	def is_member(self, x):
		for key,value in self.boundaries.items():
			foo = x.scalar(key)
			if foo <= value[0] or foo >= value[1]:
				return False
			for r in positive_roots(R):
				if x.scalar(r).denominator() == 1:
					return False
		return True

	def cover(self, x):
		R = self.root_system
		foo = fundamental_coweights(R).list()
		foo.append(sum(self.quotient_basis))
		return len([i for i in foo if self.is_member(x-i)])

	def proj_cover(self,x):
		R = self.root_system
		foo = fundamental_coweights(R).list()
		n = len(R.index_set())
		for i in R.index_set():
			foo[i-1] = foo[i-1]/(n-1) - sum(self.quotient_basis)
		return len([i for i in foo if self.is_member(x-i)])

def i_order_leq(i,n,a,b):
	if i <= a and a <= b:
		return True
	if i <= a and b <= i-1:
		return True
	if a <= b and b <= i-1:
		return True
	return False

class DecoratedPermutation:
	def __init__(self, pi, col):
		self.permutation = Permutation(pi)
		self.n = len(pi)
		self.color = col
		self.k = -1
		
	def fixed_points(self): 
		return self.permutation.fixed_points()
	
	def weak_excedance(self,i):
		wei = set()
		fp = self.fixed_points()
		for j in range(1,self.n+1):
			if (j not in fp) and i_order_leq(i,self.n,j,self.permutation(j)):
				wei.add(j)
			if j in fp and self.color[j-1] == 1:
				wei.add(j)
		self.k = len(wei)
		return wei

class GrassmannNecklace:
	def __init__(self, I):
		self.n = len(I)
		self.k = len(I[0])
		self.necklace = I
		self.i_sorted_necklace = []

	def i_sort(self):
		I = self.necklace
		n = self.n
		temp = []
		for i in range(n):
			foo = sorted([(j-i)%n for j in I[i]])
			temp.append(tuple([(j+i)%n if (j+i)%n != 0 else n for j in foo]))
		self.i_sorted_necklace = temp
		return temp

def DP_to_GN(dp):
	return GrassmannNecklace([dp.weak_excedance(i) for i in range(1,dp.n+1)])

def GN_to_DP(gn):
	n = gn.n
	I = gn.necklace
	w = []
	col = []
	for i in range(n):
		if I[(i+1)%n] == I[i]:
			w.append(i)
			if i not in I[i]:
				col.append(-1)
			else:
				col.append(1)
		else:
			w.append(I[i] - I[(i+1)%n])
			col.append(0)
	return DecoratedPermutation(w, col)

def constant_of_SR(eq):
	return eq.polynomial(ZZ).constant_coefficient()

def extract_coefficients(eq, vs):
	return [constant_of_SR(eq)] + [eq.coefficient(v) for v in vs]

def GenerateMatrix(eqsys, vs):
	A = []
	for eq in eqsys:
		RHS = eq.rhs()
		LHS = eq.lhs()
		if eq.operator() == (var('t') <= 0).operator() or eq.operator() == (var('t') < 0).operator():
			foo1 = RHS - LHS
			A.append(extract_coefficients(foo1, vs))
		if eq.operator() == (var('t') >= 0).operator() or eq.operator() == (var('t') > 0).operator():
			foo2 = LHS - RHS
			A.append(extract_coefficients(foo2, vs))
		if eq.operator() == (var('t') == 0).operator():
			foo1 = RHS - LHS
			foo2 = LHS - RHS
			A.append(extract_coefficients(foo1, vs))
			A.append(extract_coefficients(foo2, vs))
	return matrix(ZZ, A)

def GN_to_ineqs(gn):
	n = gn.n
	I = gn.i_sort()
	k = gn.k
	ineqs = []
	v = [var('x%d' % i) for i in range(1,n+1)]
	ineqs.append(sum(v) == k)
	for xi in v:
		ineqs.append(xi >= 0)
		ineqs.append(xi <= 1)
	for i in range(1,n+1):
		for j in (1..k):
			if I[i-1][j-1] - i > j-1:
				ineqs.append(sum([var('x%d' %i) for i in range(i,I[i-1][j-1])]) <= j-1)
	return ineqs, v

def GN_to_bdp(gn):
	bdp = {}
	I = gn.i_sort()
	n = gn.n
	k = gn.k
	R = RootSystem(['A', n-1])
	sr = R.root_lattice().simple_roots()
	theta = R.root_lattice().highest_root()
	for i in R.index_set():
		bdp[sr[i]] = [0,1]
	bdp[theta] = [k-1,k]
	for i in (1..n):
		for j in (2..k):
			if I[i-1][j-1] - i > j-1:
				bdp[sum([sr[k] for k in (i..I[i-1][j-1]-1)])] = [0,j-1]
	return bdp

class Positroid:
	def __init__(self, gn):
		self.necklace = gn
		self.inequalities, self.variables = GN_to_ineqs(gn)
		self.n = len(self.variables)
		self.dimension = self.n-1
		self.k = gn.k
		A = GenerateMatrix(self.inequalities, self.variables)
		self.polytope = Polyhedron(ieqs = A, backend='normaliz', base_ring=ZZ)
		self.projected_polytope = self.polytope.affine_hull_projection()
		self.root_system = RootSystem(['A', self.n-1])
		self.alcoved = AlcovedPolytope(self.root_system, GN_to_bdp(self.necklace))

	def bases(self):
		return {tuple([i+1 for i in range(self.n) if v[i] != 0]) for v in self.polytope.vertices()}
