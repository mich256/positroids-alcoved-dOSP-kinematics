# Generate all positroids on [n] from decorated permutations

def gen_dec_perm(n):
	temp = []
	#one = Permutations(n).identity()
	for pi in CyclicPermutations(range(1,n+1)):
		pi = Permutation(pi)
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