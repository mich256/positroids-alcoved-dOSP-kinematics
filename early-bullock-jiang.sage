def hDOSP(n,k):
	temp = []
	for c2 in Compositions(n,min_part=2):
		for osp in OrderedSetPartitions(n,c2):
			if 1 in osp[0]:
				for c in Compositions(k,inner=[1]*len(osp), outer=[len(i) for i in osp]):
					temp.append((osp,tuple(c)))
	return temp

def hDOSP_pp(pair):
	return ', '.join(str(set(pair[0][i]))+str(str(pair[1][i])) for i in range(len(pair[1])))

def J_Sr(J,n):
	S = []
	r = []
	m = min(i for i in range(1,n+1) if i not in J)
	i = m
	t1 = set()
	t2 = 0
	t = 0
	while i < n+1:
		if i in J:
			t2 += 1
			t = 1
		if i not in J:
			if t == 1:
				S.append(t1)
				t1 = set()
				r.append(t2)
				t2 = 0
				t = 0
		t1.add(i)
		i += 1
	t1.update(range(1,m))
	S.append(t1)
	r.append(t2)
	return (OrderedSetPartition(S),tuple(r))

def sJ_XSr(Sr,J):
	return max([sum(Sr[1][:j]) - sum(len(J.intersection(Set(Sr[0][i]))) for i in range(j)) for j in range(len(Sr[0]))])

def XSr(Sr):
	n = Sr[0].base_set_cardinality()
	k = sum(Sr[1])
	return {J: sJ_XSr(Sr,J) for J in Subsets(n,k)}

def sJ(J,n):
	temp = dict()
	k = len(J)
	IJ = set([j for j in range(1,n+1) if j in J and j-1 not in J])
	if 1 in J and n not in J:
		IJ.add(1)
	for M in Subsets(IJ):
		JJ = set()
		for j in range(1,n+1):
			if j in (J-M):
				JJ.add(j)
			if j in (J&M):
				if j != 1:
					JJ.add(j-1)
				else:
					JJ.add(n)
		temp[frozenset(JJ)] = (-1)^(len(J.intersection(M))-k-1)
	return temp

def XJ_XSr(Sr,J):
	n = Sr[0].base_set_cardinality()
	k = sum(Sr[1])
	IJ = set([j for j in range(1,n) if j in J and j+1 not in J])
	temp = 0
	if n in J and 1 not in J:
		IJ.add(n)
	for M in Subsets(IJ):
		JJ = set()
		for j in range(1,n+1):
			if j in J.difference(M):
				JJ.add(j)
			if j in J.intersection(M):
				if j != n:
					JJ.add(j+1)
				else:
					JJ.add(1)
		temp += (-1)^(len(J.intersection(M))-k-1)*sJ_XSr(Sr,JJ)
	return temp

def main(n,k):
	for Sr in hDOSP(n,k):
		d = XSr(Sr)
		print(hDOSP_pp(Sr))
		print({J: XJ_XSr(Sr,J) for J in Subsets(n,k)},'\n')

if __name__ == '__main__':
	import sys
	main(*map(int, sys.argv[1:]))
