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
			trueP = P.projected_polytope
			h1 = trueP.h_star_vector()
			h2 = P.cover_no()
			if trueP.dim() == n-1:
				if h1 != h2:
					print(P.alcoved.boundaries, P.cover_stats())
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