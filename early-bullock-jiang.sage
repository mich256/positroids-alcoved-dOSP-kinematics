from itertools import chain

def hDOSP(n,k):
	temp = []
	for c2 in Compositions(n,min_part=2):
		for osp in OrderedSetPartitions(n,c2):
			if 1 in osp[0]:
				for c in Compositions(k,inner=[1]*len(osp), outer=[len(i) for i in osp]):
					temp.append((osp,c))
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
	return (S,r)

def sJ_XSr(Sr,J):
	return max([sum(Sr[1][:j]) - sum(len(J.intersection(Set(Sr[0][i]))) for i in range(j)) for j in range(len(Sr[0]))])

def XSr(Sr):
	n = Sr[0].base_set_cardinality()
	k = sum(Sr[1])
	return {J: sJ_XSr(Sr,J) for J in Subsets(n,k)}

#def sJ(J):
