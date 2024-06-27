# Generate all positroids on [n] from decorated permutations

def gen_dec_perm(n):
	temp = []
	for pi in CyclicPermutations(n):
		pi = Permutation(pi)
		fp = pi.fixed_points()
		for binv in range(2**len(fp)):
			col = [0]*n
			binv = bin(binv)[2:]
			for i in fp:
				col[i] = binv[i]
			temp.append(DecoratedPermutation(pi, col))
	return temp

# Generate random alcoved polytopes

