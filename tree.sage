def alignments(pi):
    '''kudos to Christian Gaetz'''
    n=len(pi)
    count=0
    for i in range(1,n+1):
        for j in range(1,n+1):
            if i<j and j<pi(j) and pi(j)<pi(i):
                count+=1
            if j<pi(j) and pi(j)<pi(i) and pi(i)<i:
                count+=1
            if pi(j)<pi(i) and pi(i)<i and i<j:
                count+=1
            if pi(i)<i and i<j and j<pi(j):
                count+=1
    return count

def trees(n):
    foo = []
    for w in Derangements(n):
        w = Permutation(w)
        k = n-len(w.weak_excedences())
        if len(w.cycle_type()) == 1 and alignments(w) == (k-1)*(n-k-1):
            foo.append(w)
    return foo

def grassmann_necklace(w):
    n = len(w)
    c1 = Permutation(list(range(2,n+1))+[1])
    c2 = Permutation([n]+list(range(1,n)))
    v = Permutation(w)
    I = []
    for i in range(n):
        I.append(tuple([(v(x)+i)%n if v(x) != n-i else n for x in v.inverse().weak_excedences()]))
        v = v.left_action_product(c1).right_action_product(c2)
    return I

def gn_to_ineqs(gn):
    n = len(gn)
    I = gn
    k = len(gn[0])
    ineqs = []
    v = [var('x%d' % i) for i in range(1,n+1)]
    ineqs.append(sum(v) == k)
    for xi in v:
        ineqs.append(xi >= 0)
        ineqs.append(xi <= 1)
    for i in range(1,n+1):
        for j in (1..k):
            if (I[i-1][j-1]-i)%n > j-1:
                y = 0
                for l in range((I[i-1][j-1]-i)%n):
                    if l+i == n:
                        y += var('x%d'%n)
                    else:
                        y += var('x%d'%((l+i)%n))
                ineqs.append(y <= j-1)
    return ineqs, v

load('alcoved.sage')

def test(n):
    for w in trees(n):
        gn = grassmann_necklace(w)
        k = len(gn[0])
        ieqs, v = gn_to_ineqs(gn)
        A = matrix_from_ieqs(ieqs, v)
        P = Polyhedron(ieqs = A, backend='normaliz', base_ring=ZZ)
        d = dict((i,[(i+1)%n]) for i in range(n))
        for H in P.Hrepresentation():
            H = H.A()
            if sum(H) == n:
                continue
            l = [i for i in range(n) if H[i] != 0]
            if len(l) == 1:
                continue
            d[min(l)].append(max(l)+1)
        Graph(d).show(vertex_color='white')
        return w, factor(P.ehrhart_polynomial())
