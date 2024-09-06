load('alcoved.sage')

def subC(W):
	l = []
	for w in W:
		if cdes_WG(W,w) == 1:
			l.append(w)
	return l

# need to take s_0 into account
def genGraph(W,k):
	C = subC(W)
	d = {}
	for w in W:
		if cdes_WG(W,w.inverse()) == k:
			d[w] = []
			for s in W.simple_reflections():
				if cdes_WG(W, s*w.inverse()) == k:
					d[w].append(w*s)
	G = Graph(d)
	G.plot().show()
	return Graph(d)

def hypersimplex(W,k):
	R = RootSystem(W)
	theta = R.root_space().highest_root()
	bdp = {theta: [k-1,k]}
	for sr in R.root_space().simple_roots():
		bdp[sr] = [0,1]
	return AlcovedPolytope(R, bdp)

load('hstar.sage')