# Generate all positroids on [n] from decorated permutations
load('alcoved.sage')

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

def split_binary_string_to_int(binary_string):
    # Ensure the input is a string and contains only 0's and 1's
    if not isinstance(binary_string, str) or not all(char in '01' for char in binary_string):
        raise ValueError("Input must be a string containing only '0' and '1'.")

    # Convert the string to a list of integers
    return [int(char) for char in binary_string]

from itertools import product

def generate_binary_vectors(length):
  """Generates all binary vectors of a given length."""
  return list(product([0, 1], repeat=length))