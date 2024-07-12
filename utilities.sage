# Generate all positroids on [n] from decorated permutations

def gen_dec_perm(n):
	temp = []
	#one = Permutations(n).identity()
	for pi in Permutations(n):
		#if pi == one:
		#	continue
		fp = pi.fixed_points()
		for binv in range(2**len(fp)):
			col = [0]*n
			binv = bin(binv)[2:]
			binv = (len(fp)-len(binv))*'0' + binv
			for i in range(len(fp)):
				col[fp[i]-1] = int(binv[i])
			temp.append(DecoratedPermutation(pi, col))
	return temp

def excedences(w):
	n = len(w)
	foo = []
	for i in range(1,n+1):
		if w(i) > i:
			foo.append(i)
	return foo

def rev_exc(w):
	n = len(w)
	foo = []
	for i in range(1,n+1):
		if w(i) < i:
			foo.append(i)
	return foo

def test1(n):
	load('utilities.sage')
	load('hstar.sage')
	load('alcoved.sage')
	S = gen_dec_perm(n)
	count = 0
	for dec in S:
		G = DP_to_GN(dec)
		P = Positroid(G)
		if P.polytope.dim() > 0:
			h1 = P.polytope.h_star_vector()
			h2 = P.cover_no()
			if P.polytope.dim() == n-1:
				try:
					P1 = P.polytope.boundary_complex()
					print(P1.is_shellable(certificate=True))
				except:
					print('not simplicial \n', dec)
				#print(P1.vertices())
				if h1 != h2:
					print(False, h1, h2)
				else:
					count += 1
			#else:
				#print(dec, G.i_sorted_necklace, trueP.dim())
	return count

def test2(n):
	load('alcoved.sage')
	load('utilities.sage')
	for pi in Derangements(n):
		pi = pi.to_permutation()
		w = pi.fundamental_transformation()
		bar = set()
		for i in excedences(pi.inverse()):
			bar.add(pi.inverse()(i))
		foo = set()
		for i in w.descents():
			foo.add(w(i))
		foo = frozenset(foo)
		bar = frozenset(bar)
		if foo != bar:
			return False
		return True

def genf_des_ides(n):
	load('des_cover.sage')
	foo = 0
	R.<u,t> = PolynomialRing(ZZ,'u,t')
	#for w in CyclicPermutations(range(1,n+1)):
	#	foo += u^(cdes(w)) * t^(icdes(w))
	for w in Permutations(n):
		foo += u^(w.number_of_descents()) * t^(w.number_of_idescents())
	return foo

def A_nk(n,k):
	foo = 0
	R.<t> = PolynomialRing(ZZ,'t')
	for w in Permutations(n):
		if w.number_of_descents() == k:
			foo += t^(w.number_of_idescents())
	return foo
