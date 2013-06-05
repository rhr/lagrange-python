import sys, input, igraph, scipy, nchoosem
from pprint import pprint

class AreaGraph(igraph.Graph):
    def are_connected(self, vertices):
        return self.subgraph(vertices).is_connected()

    def subranges(self, vertices, n):
        m = len(vertices)
        assert n < m
        return [ tuple(x) for x in nchoosem.subsets(vertices, n) \
                 if self.are_connected(x) ]


## def parse_lower_triangle(s):
##     lines = s.strip().split("\n")
##     nareas = len(lines)
##     dm = scipy.ones((nareas, nareas))
##     edges = []
##     for ai, line in enumerate(lines):
##         v = line.split()[:-1]
##         for i, x in enumerate(v):
##             x = float(x)
##             dm[ai,i] = x; dm[i,ai] = x
##             if x > 0:
##                 edges.append((ai,i))
##     return dm, igraph.Graph(edges=edges)



if __name__ == "__main__":
    s = """
    -
    1 -
    1 1 -
    0 0 0 -
    0 0 0 0 -
    1 0 0 1 0 -
    1 1 0 0 0 0 -
    0 0 0 0 0 0 0 -
    0 0 0 0 0 0 0 1 -
    0 0 1 1 0 0 0 0 0 -
    0 0 0 0 0 0 0 0 0 0 -
    0 0 0 0 0 0 0 0 0 0 1 -
    0 0 0 0 0 0 0 0 0 0 1 1 -
    0 0 0 0 0 0 0 0 0 0 1 1 1 -
    1 0 0 0 0 0 0 0 0 0 0 0 0 0 -
    """

    dm, g = input.parse_lower_triangle(s)
    #print list(dm[0])
    #pprint([ e.tuple for e in g.es ])
    for i in range(len(dm)-1):
        print i, g.are_connected((0,i))

