def newton_identities(p):
	e =  p.coefficients(sparse = False)
	n = p.degree()
	e += [0]*n
	m = [n,-e[1]]
	for k in range(2,2*n-1):
		m.append(-k*e[k]-sum([e[k-i]*m[i] for i in range(1,k)]))
	return m

def hermite_matrix(p):
	m = newton_identities(p)
	n = p.degree()
	H = zero_matrix(n,n)
	for i in range(n):
		for j in range(n):
			H[i,j] = m[i+j]
	return H

def principal_submatrix(m, k):
    return m[range(k), range(k)]

def principal_minors(m):
	n = m.ncols()
	return [det(principal_submatrix(m, k)) for k in range(1,n+1)]

def real_rooted(p):
	return all([d >= 0 for d in principal_minors(hermite_matrix(p))])