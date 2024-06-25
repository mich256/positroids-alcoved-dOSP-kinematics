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

def no_hypersimplicial_dosp(n,k):
    return {key: len(value) for (key, value) in hypersimplicial_dosp(n,k).items()}

# rewrite everything in class language. object oriented!

class Family:
    def __init__(self, *args):
        if len(args) == 1:
            self.colored_sequence = args[0]
        if len(args) == 2:
            t = sorted(args[0])
            self.colored_sequence = t[args[1]:]+list(reversed(t[0:args[1]]))

    def underlying_set(self):
        return set(self.colored_sequence)
    
    def maximum(self):
        return max(self.underlying_set())
    
    def decoration(self):
        l = self.length()
        m = self.maximum()
        return l - self.colored_sequence.index(self.maximum())
    
    def anchor(self):
        return min(self.underlying_set)
    
    def sorted_list(self):
        return sorted(self.underlying_set())
    
    def lowest_red(self):
        return self.colored_sequence[0]
    
    def highest_blue(self):
        return self.sorted_list()[self.decoration()-1]

    def regular_insert(self, x):
        if x < self.highest_blue:
            return Family(self.underlying_set|{x}, self.decoration+1)
        else:
            return Family(self.underlying_set|{x}, self.decoration)

    def is_singlet(self):
        if len(self.underlying_set) == 1:
            if self.decoration != 0:
                raise Exception('singlet has no descent or ascent')
            return True
        else:
            return False

    def remove(self,x):
        self.colored_sequence.remove(x)

class FamilyRegistry:
    def __init__(self, listOfFamilies):
        self.registry = listOfFamilies

    def underlying_set(self):
        foo = [f.underlying_set for f in self.registry]
        return set().union(*foo)

    def registry_to_permutation(self):
        n = max(self.underlying_set())
        R = self.registry
        for i in range(len(R)):
            if n in R[i]:
                max_index = i
                break
        F = R[max_index]
        # Case A:
        if len(F) == 1:
            if max_index != len(R):
                raise Exception('Sorry, singlet of the max number should be at the end')
            R.pop()
            return registry_to_permutation(FamilyRegistry(R))+[n]

        # Case C:
        if len(F) > 2:
            index_of_n_in_F = F.index(n)
            friend_of_n = F[index_of_n_in_F+1]
            F = F.remove(n)
            R[max_index] = F

        if len(F) == 2:
            friend_of_n = F[1]
            if max_index == 0 or friend_of_n > min(R[max_index-1]):
                slid_right = max_index
                while slid_right < len(R)-1:
                    if friend_of_n < min(R[slid_right+1]):
                        slid_right += 1
                    else:
                        break
                if slid_right == len(R):
                    # Case B
                    R.pop(max_index)
                    R = R + Family([friend_of_n])

                else:
                    # Case E
                    R.pop(max_index)
                    R[slid_right+1] = R[slid_right+1].regular_insert(friend_of_n)

            else:
                # Case D
                slid_left = max_index
                while slid_left > 0:
                    if friend_of_n < min(R[slid_left-1]):
                        slid_left -= 1
                    else:
                        break
                N = R[slid_left+1].regular_insert(friend_of_n)
                R[max_index] = N
                R.pop(slid_left+1)

        pi = registry_to_permutation(FamilyRegistry(R))
        index_of_friend_in_pi = pi.index(friend_of_n)
        return pi.insert(index_of_friend_in_pi-1, n)

def permutation_to_registry(w):
    n = max(w)
    max_index = w.index(n)
    if max_index == len(w):
        # Case A
        w.remove(n)
        R = permutation_to_registry(w)
        return FamilyRegistry(R+[Family([n])])
    m = w[max_index+1] # friend of n
    w.remove(n)
    R = permutation_to_registry(w).registry
    for i in len(R):
        if m in R[i]:
            index_of_m = i
            break
    F = R[index_of_m]
    # Case B and C
    if m == F.highest_blue():
        F = F.regular_insert(n)
        R[index_of_m] = F
    else:
        N = Family([n,m])
        # Case D
        if m == F.anchor():
            R[index_of_m] = N
            slid_left = index_of_m
            while slid_left > 0:
                if m < min(R[slid_left-1]):
                    slid_left -= 1
                else:
                    break
            F.remove(m)
            R.insert(slid_left, F)
        # Case E
        else:
            F.remove(m)
            R[index_of_m] = F
            slid_left = index_of_m
            while slid_left > 0:
                if m < min(R[slid_left-1]):
                    slid_left -= 1
                else:
                    break
            R.insert(slid_left, N)

    return FamilyRegistry(R)

