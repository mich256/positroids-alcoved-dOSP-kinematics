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

def inverse_circular_descents(w):
	return circular_descents(w.inverse())

def cdes(w):
	w = Permutation(w)
	return len(circular_descents(w))

def icdes(w):
	w = Permutation(w)
	return cdes(w.inverse())

def cover(w):
	w = Permutation(w)
	n = len(w)
	foo = [i for i in range(1,n) if w(i+1)>w(i)+1]
	temp = len(foo)
	if w(1) == 1:
		return temp
	else:
		return temp+1

def circular_cover(w):
	w = Permutation(w)
	n = len(w)
	foo = []
	for i in range(1,n):
		if w(i) == 1 and w(i+1) == n:
			# this will result in lower cdes in general
			continue
		elif w(i)+1 < w(i+1):
			foo.append(i)
	return len(foo)

def perm_icdes(n,k):
	return [Permutation(w) for w in CyclicPermutations(range(1,n+1)) if icdes(w) == k]

def perm_ides(n,k):
	return [w for w in Permutations(n) if w.number_of_idescents() == k]

def perm_cover(n,k):
	temp = {}
	for w in perm_ides(n,k):
		c = cover(w)
		temp.setdefault(c, [])
		temp[c].append(w)
	return temp

def perm_c_cover(n,k):
	temp = {}
	for w in perm_icdes(n,k):
		c = circular_cover(w)
		temp.setdefault(c, [])
		temp[c].append(w)
	return temp

def no_perm_c_cover(n,k):
	return {key: len(value) for (key, value) in perm_c_cover(n,k).items()}

def perm_c_cover_polynomial(n,k):
	R.<t> = PolynomialRing(QQ)
	temp = no_perm_c_cover(n,k)
	return sum([value*t^key for (key, value) in temp.items()])

def no_perm_cover(n,k):
	return {key: len(value) for (key, value) in perm_cover(n,k).items()}

def perm_cover_polynomial(n,k):
	R.<t> = PolynomialRing(QQ)
	temp = no_perm_cover(n,k)
	return sum([value*t^key for (key, value) in temp.items()])

load('hstar.sage')

def test_c_cover(n,k):
	R.<t> = PolynomialRing(QQ)
	return perm_c_cover_polynomial(n,k) == h_star_polynomial_of_hypersimplex(n,k)

def test_cover(n,k):
	R.<t> = PolynomialRing(QQ)
	return perm_cover_polynomial(n-1,k-1) + (1-t)*h_star_polynomial_of_hypersimplex(n-1,k-1) \
	== h_star_polynomial_of_hypersimplex(n,k)