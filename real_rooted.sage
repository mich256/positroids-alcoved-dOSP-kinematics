def newton_identities(p):
	e =  p.coefficients(sparse = False)
	n = len(e)-1
	m = [n,-e[1]]
	for k in range(2,n+1):
		m.append(-k*e[k]-sum([e[k-i]*m[i] for i in range(1,k)]))
	return m

def Hmatrix(p):
	m = newton_identities(p)
	n = len(m)-1
	H = zero_matrix(n,n)
	for i in range(n):
		for j in range(n):
			if i+j < n+1:
				H[i,j] = m[i+j]
	return H