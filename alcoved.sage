load('des_cover.sage')
load('decorated_osp.sage')

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

def coweight_to_permutation(R, coweight):
	n = len(R.index_set())
	w = [1]+[0]*(n-1)
	sr = R.root_lattice().simple_roots()
	bar = 0
	count = 2
	for i in sr:
		foo = coweight.scalar(i)
		w[(bar + i)%n] = count
		bar = bar + i
		count += 1
	return Permutation(w)

def permutation_to_coweight(u):
	n = len(u)
	R = RootSystem(['A', n-1])
	lcheck = fundamental_coweights(R)
	w = u.inverse()
	foo = 0
	for i in (1..n-1):
		foo += Rational(((w(i+1)-w(i))%n)/n) * lcheck[i]
	return foo

def R_to_W(R, alpha):
	W = WeylGroup(R.cartan_type())
	A = W.domain().simple_roots()
	return sum([alpha.coefficient(i) * A[i] for i in R.index_set()])

def inv(R, w, alpha):
	W = WeylGroup(R.cartan_type())
	A = W.domain().simple_roots()
	if w.action(alpha).is_positive_root():
		return 0
	if (-w.action(alpha)).is_positive_root():
		return 1
	else:
		raise Exception('neither positive nor negative')

def cdes_wrt_root(R, w, alpha):
	W = WeylGroup(R.cartan_type())
	A = W.domain().simple_roots()
	foo = 0
	for i in R.index_set():
		foo += alpha.coefficient(i) * inv(R, w, A[i])
	foo += inv(R, w, -R_to_W(R,alpha))
	return foo

def cdes_WG(R, w):
	theta = R.root_lattice().highest_root()
	return cdes_wrt_root(R, w, theta)

def permutation_to_typeA_WeylGroup(w):
	n = len(w)
	W = WeylGroup(['A',n-1])
	s = W.simple_reflections()
	foo = w.reduced_word()
	bar = W.one()
	for i in foo:
		bar *= s[i]
	return bar

def des_to_bin(d,n):
	foo = [0]*n
	for i in d:
		foo[i-1] = 1
	return ''.join(map(str,foo))

def circuits(w):
	n = len(w)
	cycle = list(range(2,n+1))+[1]
	cycle = Permutation(cycle)
	foo = []
	for i in range(n):
		#foo.append(des_to_bin(circular_descents(w.inverse()),n))
		foo.append(frozenset(circular_descents(w.inverse())))
		w = w.left_action_product(cycle)
	return foo

def circuits_of_WG(w):
	#TODO
	return

def BDP_to_ieqs(R, BDP):
	sr = R.ambient_space().simple_roots()
	ieqs = []
	n = R.ambient_space().dimension()
	v = [var('x%d' % i) for i in (1..n-1)]
	for key, value in BDP.items():
		foo = 0
		for i in R.index_set():
			foo += key.coefficient(i) * sr[i]
		ieqs.append(sum([foo[i-1] * var('x%d' %i) for i in (1..(n-1))]) >= value[0])
		ieqs.append(sum([foo[i-1] * var('x%d' %i) for i in (1..(n-1))]) <= value[1])
	return ieqs, v


class AlcovedPolytope:
	def __init__(self, R, BoundaryParameters):
		self.root_system = R
		self.boundaries = BoundaryParameters
		# the coweights are indexed from 1, not from 0
		self.quotient_basis = [fundamental_coweights(R)[i]/R.cartan_type().dual_coxeter_number() for i in R.index_set()]
		ieqs, v = BDP_to_ieqs(R, BoundaryParameters)
		A = GenerateMatrix(ieqs, v)
		self.polytope = Polyhedron(ieqs = A, backend='normaliz', base_ring=QQ)

	def is_member(self, x):
		R = self.root_system
		for key,value in self.boundaries.items():
			foo = x.scalar(key)
			if foo <= value[0] or foo >= value[1]:
				return False
			for r in positive_roots(R):
				if x.scalar(r).denominator() == 1:
					return False
		return True


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

	def __repr__(self):
		fp = self.fixed_points()
		foo = ''
		for i in range(1,self.n+1):
			if i in fp:
				if self.color[i-1] == 1:
					foo += str(self.permutation(i))+'~'
				else:
					foo += str(self.permutation(i))+'_'
			else:
				foo += str(self.permutation(i))
		return foo

	def __str__(self):
		return self.__repr__()
	
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

	def __repr__(self):
		return str(self.necklace)
	
	def __str__(self):
		return str(self.necklace)

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
	return GrassmannNecklace([frozenset(dp.weak_excedance(i)) for i in range(1,dp.n+1)])

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
			if (I[i-1][j-1]-i)%n > j-1:
				y = 0
				for l in range((I[i-1][j-1]-i)%n):
					if l+i == n:
						y += var('x%d'%n)
					else:
						y += var('x%d'%((l+i)%n))
				ineqs.append(y <= j-1)
	return ineqs, v

def GN_to_bdp(gn):
	bdp = {}
	I = gn.i_sort()
	n = gn.n
	k = gn.k
	R = RootSystem(['A',n-1])
	sr = R.root_lattice().simple_roots()
	theta = R.root_lattice().highest_root()
	for i in R.index_set():
		bdp[sr[i]] = [0,1]
	bdp[theta] = [k-1,k]
	for i in (1..n):
		#for j in (2..k):
		for j in (1..k):
			if (I[i-1][j-1] - i)%n > j-1:
				y = 0
				for l in range((I[i-1][j-1]-i)%n):
					if l+i == n:
						y -= theta
					else:
						y += sr[(l+i)%n]
				if y in positive_roots(R):
					bdp[y] = [0,j-1]
				else:
					foo = -y
					if foo in bdp.keys():
						value = bdp[foo]
						bdp[foo] = [max(k-j+1,value[0]), min(len(foo.coefficients()),value[1])]
					else:
						bdp[foo] = [k-j+1, len(foo.coefficients())]
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
		if self.polytope.dim() >= 0:
			self.projected_polytope = self.polytope.affine_hull_projection()
		self.root_system = RootSystem(['A', self.n-1])
		self.alcoved = AlcovedPolytope(self.root_system, GN_to_bdp(self.necklace))

	def is_member(self,w):
		x = permutation_to_coweight(w)
		if self.alcoved.is_member(x):
			return True
		return False

	def all_members(self):
		foo = []
		for w in perm_icdes(self.n,self.k):
			if self.is_member(w):
				foo.append(w)
		return foo

	def cover(self, s, w):
		n = self.n
		#P = Permutations(self.n)
		#s = P.simple_reflections()
		foo = 0
		if self.is_member(w):
			for i in w.descents():
				if i == n-1 and w(i) == n:
					continue
				else:
					u = w.left_action_product(s[i])
					a = w(i)
					b = w(i+1)
					if self.is_member(u) and icdes_ij(u,a,b) < icdes_ij(w,a,b):
						foo += 1
		return foo

	def cover_stats(self):
		P = Permutations(self.n)
		s = P.simple_reflections()
		foo = {}
		for i in range(self.n+1):
			foo[i] = []
		for w in self.all_members():
			foo[self.cover(s,w)].append(w)
		return foo

	def cover_no(self):
		foo = self.cover_stats()
		bar = []
		for i in range(self.n+1):
			if len(foo[i]) == 0:
				break
			else:
				bar.append(len(foo[i]))
		return bar

	def bases(self):
		return {tuple([i+1 for i in range(self.n) if v[i] != 0]) for v in self.polytope.vertices()}
