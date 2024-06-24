class AlcovedPolytope:
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
		if self.is_member(coweight) == False:
			raise Exception('coweight is not in the alcoved polytope')
		foo = self.fundamental_coweights.list().append(self.rho)
		return sum([i for i in foo if self.is_member(coweight-i)])

def i_order_eq(i,n,a,b):
	if i <= a and a <= b:
		return True
	if i <= a and b <= i-1:
		return True
	if a <= b and b <= i-1:
		return True
	return False

class DecoratedPermutation:
	def __init__(self, pi, col):
		self.permutation = pi
		self.n = len(pi)
		self.k = None
		self.color = col
		self.fixed_points = pi.fixed_points()

	def weak_excedance(self,i):
		wei = set()
		for j in range(n):
			if j not in self.fixed_points and i_order_eq(i,n,j,self.permutation(j)):
				wei.add(j)
			if j in self.fixed_points and self.color[j] == 1:
				wei.add(j)
		self.k = len(wei)
		return wei

class GrassmannNecklace:
	def __init__(self, I):
		self.n = len(I)
		self.k = len(I[0])
		self.necklace = I
		self.sorted_necklace = None

def DP_tp_GN(dp):
	return [dp.weak_excedance(i) for i in range(dp.n)]

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

def i_sort(self):
	I = self.necklace
	n = self.n
	temp = []
	for i in range(n):
		foo = sorted([(j-i)%n for j in I[i]])
		temp.append(((j+i)%n for j in foo))
	self.sorted_necklace = temp
	return temp

def GenerateMatrix(equsys, vs):
	A = []
	for eq in equsys:
		if eq.operator() == (x <= 0).operator() or eq.operator() == (x < 0).operator():
			A.append((equ.rhs())+(-equ.lhs().coefficient(v) for v in vs))
		if eq.operator() == (x >= 0).operator() or eq.operator() == (x > 0).operator():
			A.append((-equ.rhs())+(equ.lhs().coefficient(v) for v in vs))
		if eq.operator() == (x == 0).operator():
			A.append((-equ.rhs())+(equ.lhs().coefficient(v) for v in vs))
			A.append((equ.rhs())+(-equ.lhs().coefficient(v) for v in vs))
	return A

def GN_to_ineqs(gn):
	n = gn.n
	gn.i_sort()
	I = gn.sorted_necklace
	k = gn.k
	ineqs = []
	v = [var('x%d' % i) for i in range(n)]
	ineqs.append(sum([var('x%d' % i) for i in range(n-1)] <= k))
	ineqs.append(sum([var('x%d' % i) for i in range(n-1)] >= k-1))
	for i in range(n-1):
		ineqs.append(var('xi') >= 0)
	for i in range(n):
		for j in range(k):
			ineqs.append(sum([var('x%d' %i) for i in range(i,I[i][j])] <= j-1))
	return ineqs, v

class PositroidPolytope:
	def __init__(self, gn):
		self.necklace = gn
		self.inequalities, self.variables = GN_to_ineqs(gn)
		self.polytope = Polyhedron(ineqs = [GenerateMatrix(self.inequalities, self.variables)])
