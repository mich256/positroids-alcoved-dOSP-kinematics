def circular_descents(w):
	n = len(w)
	return [i for i in range(n) if w[i] > w[(i+1)%n]]

def cdes(w):
	return len(circular_descents(w))

def des(w):
	n = len(w)
	return len([i for i in range(n-1) if w[i] > w[i+1]])

def cover(w):
	n = len(w)
	temp = sum([int(w[(i+1)%n]>w[i]+1) for i in range(n)])
	if w(1) == 1:
		return temp
	else:
		return temp+1

def perm_cdes(n,k):
	return [w for w in CyclicPermutations(range(1,n+1)) if cdes(w) == k]

def perm_des(n,k):
	return [w for w in Permutations(n) if w.number_of_descents() == k]

def perm_cover(n,k):
	temp = {}
	for w in perm_des(n,k):
		c = cover(w)
		temp.setdefault(c, [])
		temp[c].append(w)
	return temp

def no_perm_cover(n,k):
	return {key: len(value) for (key, value) in perm_cover(n,k).items()}