def circular_descents(w):
	w = Permutation(w)
	n = len(w)
	foo = []
	for i in range(1,n):
		if w(i) > w(i+1):
			foo.append(i)
	if w(n) > w(1):
		foo.append(n)
	return foo

def cdes_ij(w, i, j):
	if i == j:
		return 0
	foo = 0
	if i < j:
		for k in range(i,j):
			if w(k) > w(k+1):
				foo += 1
		if w(j) > w(i):
			foo += 1
		return foo
	if i > j:
		n = len(w)
		for k in range(i,n):
			if w(k) > w(k+1):
				foo += 1
		if w(n) > w(1):
			foo += 1
		if j > 1:
			for k in range(1,j):
				if w(k) > w(k+1):
					foo += 1
		if w(j) > w(i):
			foo += 1
		return foo

def icdes_ij(w, i, j):
	u = w.inverse()
	return cdes_ij(u,i,j)

def ascents(w):
	w = Permutation(w)
	n = len(w)
	foo = []
	for i in range(1,n):
		if w(i) < w(i+1):
			foo.append(i)
	return foo

def inverse_circular_descents(w):
	return circular_descents(w.inverse())

def cdes(w):
	w = Permutation(w)
	return len(circular_descents(w))

def icdes(w):
	w = Permutation(w)
	u = w.inverse()
	return cdes(u)

def cover(w):
	w = Permutation(w)
	n = len(w)
	foo = [i for i in range(1,n-1) if w(i+1)>w(i)+1]
	temp = len(foo)
	if w(1) == 1:
		return temp
	else:
		return temp+1

def ncov(w):
	w = Permutation(w).inverse()
	n = len(w)
	foo = [i for i in range(1,n) if w(i+1)>w(i)+1]
	return len(foo)

def perm_ncov(n,k):
	temp = {}
	for w in perm_cdes(n,k):
		c = ncov(w)
		temp.setdefault(c, [])
		temp[c].append(w)
	return temp

def perm_icdes(n,k):
	cycle = list(range(2,n+1))
	cycle.append(1)
	cyclc = Permutation(cycle)
	return [Permutation(w).left_action_product(cycle) for w in CyclicPermutations(range(1,n+1)) if icdes(w) == k]

def perm_cdes(n,k):
	cycle = list(range(2,n+1))
	cycle.append(1)
	cyclc = Permutation(cycle)
	return [Permutation(w).left_action_product(cycle) for w in CyclicPermutations(range(1,n+1)) if cdes(w) == k]

def perm_ides(n,k):
	return [Permutation(w) for w in Permutations(n) if w.number_of_idescents() == k]

def perm_cover(n,k):
	temp = {}
	for w in perm_ides(n-1,k-1):
		c = cover(w)
		temp.setdefault(c, [])
		temp[c].append(w.inverse())
	return temp

def no_perm_cover(n,k):
	return {key: len(value) for (key, value) in perm_cover(n,k).items()}

def perm_cover_polynomial(n,k):
	R.<t> = PolynomialRing(QQ)
	temp = no_perm_cover(n,k)
	return sum([value*t^key for (key, value) in temp.items()])