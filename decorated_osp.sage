def winding(osp, decoration):
    n = osp.base_set_cardinality()
    m = len(osp)
    w = osp.to_packed_word()
    wind = []
    k = decoration.size()
    for i in range(n):
        if w[i] == w[(i+1)%n]:
            wind.append(0)
        if w[i] < w[(i+1)%n]:
            wind.append(sum([decoration[j] for j in range(w[i], w[(i+1)%n])]))
        if w[i] > w[(i+1)%n]:
            l = list(range(w[i],m))+list(range(w[(i+1)%n]))
            wind.append(sum([decoration[j] for j in l]))
    return wind, sum(wind) // k

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
                temp.setdefault(winding(osp,c)[1], [])
                temp[winding(osp,c)[1]].append([osp, c])
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
                temp.setdefault(winding(osp,c)[1], [])
                temp[winding(osp,c)[1]].append([osp, c])
    return temp

def family_registry(n,k):
    foo = {}
    bar = hypersimplicial_dosp(n,k)
    for winding_no in bar.keys():
        foo.setdefault(winding_no, [])
        for hdosp, decoration in bar[winding_no]:
            temp = []
            m = len(hdosp)
            for i in range(m):
                family = sorted(hdosp[i])
                l = len(family)
                family = family[l-decoration[i]:]+list(reversed(family[0:l-decoration[i]]))
                temp.append(family)
            foo[winding_no].append(temp)
    return foo

def no_hypersimplicial_dosp(n,k):
    return {key: len(value) for (key, value) in hypersimplicial_dosp(n,k).items()}

def registry_to_permutation(R):
    n = 0
    max_index = 0
    m = len(R)
    for i in range(m):
        family = R[i]
        if max(family) > n:
            n = max(family)
            max_index = i
    # Case A:
    if len(R[max_index]) == 1:
        R.pop()
        return registry_to_permutation(R)+[n]
    # Case C:
    if len(R[max_index]) > 2:
        index_of_n_in_F = R[max_index].index(n)
        friend_of_n = R[index_of_n_in_F+1]
        F = R[max_index].remove(n)
        R[max_index] = F
        pi = registry_to_permutation(R)
        index_of_friend_in_pi = pi.index(friend_of_n)
        return pi.insert(index_of_friend_in_pi-1, n)
    # Case B:
    TODO