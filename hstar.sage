def GenerateMatrix(equsys, vars):
    A=matrix([[equ.rhs()]+[equ.lhs().coefficient(v) for v in vars] for equ in equsys])
    return A

def h_star_vector(P):
    d = P.dim()
    E = P.ehrhart_polynomial()
    return [sum(E(t)*binomial(d+1,j-t)*(-1)^(j-t) for t in (0..j)) for j in (0..d)]

def hypercube(d):
    return Polyhedron(vertices = list(product([0,1],repeat=d)))

def Hmatrix(P):
    pH = p.Hrepresentation()
    return (pH.b().transpose() + pH.A().transpose()).transpose()

def h_star_of_hypersimplex(n,k):
    return h_star_vector(polytopes.hypersimplex(n,k))

def n_choose_2(n):
    return Subsets(list(range(n), 2))

def winding(osp, decoration):
    n = osp.base_set_cardinality()
    m = len(osp)
    w = osp.to_packed_word()
    k = decoration.size()
    return sum([(w[i+1]-w[i])%m for i in range(n-1)])//k

def CyclicOrderedSetPartitions(n):
    temp = []
    for c in OrderedSetPartitions(n):
        if 1 in c[0]:
            temp.append(c)
    return temp

def decorated_osp(n,k):
    temp = {}
    for osp in CyclicOrderedSetPartitions(n):
        m = len(osp)
        for c in Compositions(k):
            if len(c) == m:
                temp.setdefault(winding(osp,c), [])
                temp[winding(osp,c)].append([osp, c])
    return temp

def no_decorated_osp(n,k):
    return {key: len(value) for (key, value) in decorated_osp(n,k).items()}

def hypersimplicial_dosp(n,k):
    temp = {}
    for osp in CyclicOrderedSetPartitions(n):
        m = len(osp)
        for c in Compositions(k):
            if len(c) != m:
                continue
            hypersimplicial = True
            for i in range(m):
                if c[i] > len(osp[i])-1:
                    hypersimplicial = False
                    break
            if hypersimplicial:
                temp.setdefault(winding(osp,c), [])
                temp[winding(osp,c)].append([osp, c])
    return temp

def no_hypersimplicial_dosp(n,k):
    return {key: len(value) for (key, value) in hypersimplicial_dosp(n,k).items()}

