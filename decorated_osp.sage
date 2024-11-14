def CyclicOrderedSetPartitions(n):
    temp = []
    for c in OrderedSetPartitions(n):
        m = len(c)
        if 1 in c[m-1]:
            temp.append(c)
    return temp

class Block:
    def __init__(self, s, d):
        self.underlying_set = s
        self.decoration = d

    def __repr__(self):
        return str(self.underlying_set) + '_' + str(self.decoration)

    def __str__(self):
        return str(self.underlying_set) + '_' + str(self.decoration)


class CyclicDecoratedOrderedSetPartition:
    def __init__(self, *args):
        if len(args) == 1:
            self.blocks = args[0]
            self.osp = []
            self.decoration = []
            for i in range(len(self.blocks)):
                self.osp.append(self.blocks[i].underlying_set)
                self.decoration.append(self.blocks[i].decoration)
            self.osp = OrderedSetPartition(self.osp)

        if len(args) == 2:
            if len(args[0]) != len(args[1]):
                raise Exception('lengths shall be equal')
            else:
                self.osp = args[0]
                self.decoration = args[1]
                self.blocks = []
                for i in range(len(args[0])):
                    self.blocks.append(Block(args[0][i], args[1][i]))

    def descents(self):
        n = self.osp.base_set_cardinality()
        w = self.osp.to_packed_word()
        foo = []
        for i in range(n):
            if w[i] > w[(i+1)%n]:
                foo.append(i+1)
        return foo

    def number_of_descents(self):
        return len(self.descents())

    def winding(self):
        n = self.osp.base_set_cardinality()
        m = len(self.osp)
        w = self.osp.to_packed_word()
        wind = []
        k = self.decoration.size()
        for i in range(n):
            if w[i] == w[(i+1)%n]:
                wind.append(0)
            if w[i] < w[(i+1)%n]:
                wind.append(sum([self.decoration[j] for j in range(w[i], w[(i+1)%n])]))
            if w[i] > w[(i+1)%n]:
                l = list(range(w[i],m))+list(range(w[(i+1)%n]))
                wind.append(sum([self.decoration[j] for j in l]))
        return wind, sum(wind) // k

    def __repr__(self):
        return str(self.blocks)

    def __str__(self):
        return str(self.blocks)
        

def decorated_osp(n,k):
    temp = {}
    for osp in CyclicOrderedSetPartitions(n):
        m = len(osp)
        for c in Compositions(k):
            if len(c) == m:
                foo = CyclicDecoratedOrderedSetPartition(osp,c)
                bar = foo.winding()[1]
                temp.setdefault(bar, [])
                temp[bar].append(foo)
    return temp

def no_decorated_osp(n,k):
    return {key: len(value) for (key, value) in decorated_osp(n,k).items()}

# rewrite everything in class language. object oriented!

class Family:
    def __init__(self, *args):
        if len(args) == 1:
            self.colored_sequence = list(args[0])
        if len(args) == 2:
            t = sorted(args[0])
            self.colored_sequence = t[args[1]:]+list(reversed(t[0:args[1]]))

    def __str__(self):
        return str(self.colored_sequence)

    def underlying_set(self):
        return set(self.colored_sequence)
    
    def max(self):
        return max(self.underlying_set())
    
    def decoration(self):
        l = len(self.colored_sequence)-1
        if l == 0:
            return 1
        m = self.max()
        return l - self.colored_sequence.index(m)

    def __repr__(self):
        return str(self.underlying_set()) + '_' + str(self.decoration())
    
    def anchor(self):
        return min(self.underlying_set())
    
    def sorted(self):
        return sorted(self.underlying_set())
    
    def lowest_red(self):
        return self.colored_sequence[0]
    
    def highest_blue(self):
        return self.sorted()[self.decoration()-1]

    def regular_insert(self, x):
        if len(self.colored_sequence) == 1:
            raise Exception('regular insert only into nonsinglet families')
        if x < self.highest_blue():
            return Family(self.underlying_set()|{x}, self.decoration()+1)
        else:
            return Family(self.underlying_set()|{x}, self.decoration())

    def insert_with_friend(self, x, m):
        if m not in self.colored_sequence:
            raise Exception('your friend is not here')
        if x < max(self.underlying_set()):
            raise Exception('insert in natural order')
        i = self.colored_sequence.index(m)
        foo = self.colored_sequence
        foo.insert(i,x)
        return Family(foo)

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
    def __init__(self, *args):
        if len(args) == 1:
            self.registry = args[0]
            self.osp = []
            self.decoration = []
            for i in range(len(self.registry)):
                self.osp.append(self.registry[i].underlying_set())
                self.decoration.append(self.registry[i].decoration())
            self.decoration = Composition(self.decoration)
            self.osp = OrderedSetPartition(self.osp)
        if len(args) == 2:
            if len(args[0]) != len(args[1]):
                raise Exception('lengths shall be equal')
            else:
                self.osp = args[0]
                self.decoration = Composition(args[1])
                self.registry = []
                for i in range(len(args[0])):
                    self.registry.append(Family(args[0][i], args[1][i]))

    def __repr__(self):
        return str(self.registry)

    def __str__(self):
        return str(self.registry)

    def underlying_set(self):
        foo = [f.underlying_set() for f in self.registry]
        return set().union(*foo)

    def descents(self):
        n = self.osp.base_set_cardinality()
        w = self.osp.to_packed_word()
        foo = []
        for i in range(n):
            if w[i] > w[(i+1)%n]:
                foo.append(i+1)
        return foo

    def number_of_descents(self):
        return len(self.descents())

    def winding(self):
        n = self.osp.base_set_cardinality()
        m = len(self.osp)
        w = self.osp.to_packed_word()
        wind = []
        k = self.decoration.size()
        for i in range(n):
            if w[i] == w[(i+1)%n]:
                wind.append(0)
            if w[i] < w[(i+1)%n]:
                wind.append(sum([self.decoration[j] for j in range(w[i], w[(i+1)%n])]))
            if w[i] > w[(i+1)%n]:
                l = list(range(w[i],m))+list(range(w[(i+1)%n]))
                wind.append(sum([self.decoration[j] for j in l]))
        return wind, sum(wind) // k

    def to_perm(self):
        temp = []
        for f in self.registry:
            temp += f.colored_sequence
        return temp

def permutation_to_registry(w):
    if len(w) == 1:
        return FamilyRegistry([Family(w)])
    n = max(w)
    max_index = w.index(n)
    if max_index == len(w)-1:
        # Case A
        w = list(w)
        w.remove(n)
        w = Permutation(w)
        R = permutation_to_registry(w)
        return FamilyRegistry(R.registry+[Family([n])])
    m = w[max_index+1] # friend of n
    w = list(w)
    w.remove(n)
    w = Permutation(w)
    R = permutation_to_registry(w).registry
    for i in range(len(R)):
        if m in R[i].underlying_set():
            index_of_m = i
            break
    F = R[index_of_m]
    # Case B and C
    if m == F.highest_blue() or m == F.colored_sequence[0] > F.colored_sequence[1]:
        F = F.insert_with_friend(n,m)
        R[index_of_m] = F
    else:
        N = Family([n,m])
        # Case D
        if m == F.anchor():
            R[index_of_m] = N
            slid_left = index_of_m
            while slid_left > 0:
                if m < min(R[slid_left-1].colored_sequence):
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
                if m < min(R[slid_left-1].colored_sequence):
                    slid_left -= 1
                else:
                    break
            R.insert(slid_left, N)

    return FamilyRegistry(R)

def registry_to_permutation(R):
    n = max(R.underlying_set())
    R = R.registry
    if len(R) == 0 or len(R) == 1:
        return R[0].colored_sequence
    for i in range(len(R)):
        if len(R[i].colored_sequence) == 1:
            if i != len(R) - 1:
                return registry_to_permutation(FamilyRegistry(R[:i]))+R[i].colored_sequence+registry_to_permutation(FamilyRegistry(R[i+1:]))
            else:
                return registry_to_permutation(FamilyRegistry(R[:i]))+R[i].colored_sequence
    for i in range(len(R)):
        if n in R[i].underlying_set():
            max_index = i
            break
    F = R[max_index]
    # Case A:
    if len(F.colored_sequence) == 1:
        if max_index != len(R)-1:
            raise Exception('Sorry, singlet of the max number should be at the end')
        R.pop()
        return registry_to_permutation(FamilyRegistry(R))+[n]
    # Case C:
    elif len(F.colored_sequence) > 2:
        index_of_n_in_F = F.colored_sequence.index(n)
        friend_of_n = F.colored_sequence[index_of_n_in_F+1]
        F.remove(n)
        R[max_index] = F
    elif len(F.colored_sequence) == 2:
        friend_of_n = F.colored_sequence[1]
        if max_index == 0 or friend_of_n > min(R[max_index-1].underlying_set()):
            slid_right = max_index
            while slid_right < len(R)-1:
                if friend_of_n < min(R[slid_right+1].underlying_set()):
                    slid_right += 1
                else:
                    break
            if slid_right == len(R)-1:
                # Case B
                R.pop(max_index)
                R.append(Family([friend_of_n]))
            else:
                # Case E
                R[slid_right+1] = R[slid_right+1].regular_insert(friend_of_n)
                R.pop(max_index)
        else:
            # Case D
            slid_left = max_index
            while slid_left > 0:
                if friend_of_n < min(R[slid_left-1].underlying_set()):
                    slid_left -= 1
                else:
                    break
            N = R[slid_left].regular_insert(friend_of_n)
            R[max_index] = N
            R.pop(slid_left)
    pi = registry_to_permutation(FamilyRegistry(R))
    index_of_friend_in_pi = pi.index(friend_of_n)
    pi.insert(index_of_friend_in_pi, n)
    return pi

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
                foo = FamilyRegistry(osp,c)
                bar = foo.winding()[1]
                #foo = Permutation(registry_to_permutation(foo)).inverse()
                temp.setdefault(bar, [])
                temp[bar].append(foo)
    return temp

def no_hypersimplicial_dosp(n,k):
    return {key: len(value) for (key, value) in hypersimplicial_dosp(n,k).items()}


