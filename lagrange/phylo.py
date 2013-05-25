#import sets
import tree
PREORDER = 0; POSTORDER = 1
BRANCHLENGTH = 0; INTERNODES = 1

class Node:
    def __init__(self):
        self.data = {}
        self.isroot = False
        self.istip = False
        self.label = None
        self.length = None
        self.parent = None
        self.children = []
        self.nchildren = 0
        self.excluded_dists = []

    def assign_fossil(self, age, area_indices):
        t = self.age
        for seg in self.segments:
            t += seg.duration
            if t - age < 0.0001:
                # age is at beginning of existing segment
                seg.fossils = area_indices
                return
        t = self.age
        for i, seg in enumerate(self.segments):
            t0 = t + seg.duration
            if t < age < t0:
                # split this segment
                break
            t = t0
        d1 = age - t
        d2 = t0 - age
        seg1 = tree.BranchSegment(d1, seg.period, seg.model)
        seg2 = tree.BranchSegment(d2, seg.period, seg.startdist)
        self.segments = self.segments[:i] + [seg1, seg2] + self.segments[i+1:]

    def labelset_nodemap(self, d=None):
        if d is None:
            d = {}
        if not self.istip:
            s = set()
            for child in self.children:
                if child.istip:
                    s.add(child.label)
                else:
                    child.labelset_nodemap(d)
                    s = s | d[child]
            d[self] = s
            d[frozenset(s)] = self
        return d

    def mrca(self, labels):
        s = set(labels)
        tipsets = self.labelset_nodemap()
        for n in [ x for x in self.iternodes(POSTORDER) if not x.istip ]:
            if s.issubset(tipsets[n]):
                return n

    def order_subtrees_by_size(self, n2s=None, recurse=False, reverse=False):
        if n2s is None:
            n2s = node2size(self)
        if not self.istip:
            v = [ (n2s[c], c.label, c) for c in self.children ]
            v.sort()
            if reverse:
                v.reverse()
            self.children = [ x[-1] for x in v ]
            if recurse:
                for c in self.children:
                    c.order_subtrees_by_size(n2s, recurse=True, reverse=reverse)

    def add_child(self, child):
        assert child not in self.children
        self.children.append(child)
        child.parent = self
        self.nchildren += 1

    def remove_child(self, child):
        assert child in self.children
        self.children.remove(child)
        child.parent = None
        self.nchildren -= 1

##     def leaves(self, v=None):
##         if v is None:
##             v = []
##         if not self.children:
##             v.append(self)
##         else:
##             for child in self.children:
##                 child.leaves(v)
##         return v

    def leaves(self):
        return [ n for n in self.iternodes() if n.istip ]

    def iternodes(self, order=PREORDER, v=None):
        """
        returns a list of nodes descendant from self - including self
        """
        if order == PREORDER:
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER, v=None):
        """
        returns a list of nodes descendant from self - not including self!
        """
        if v is None:
            v = []
        assert order in (PREORDER, POSTORDER)
        for child in self.children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child.children:
                child.descendants(order, v)
        return v

    def find_descendant(self, label):
        if label == self.label:
            return self
        else:
            for child in self.children:
                n = child.find_descendant(label)
                if n:
                    return n
        return None

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p

    def excise(self):
        p = self.prune()
        if p.parent and len(p.children) == 1:
            c = p.children[0]
            if p.length:
                c.length += p.length
            c.prune()
            p.prune().add_child(c)
            return True
        return False

    def graft(self, node):
        parent = self.parent
        parent.remove_child(self)
        n = Node()
        n.add_child(self)
        n.add_child(node)
        parent.add_child(n)

    def leaf_distances(self, store=None, measure=BRANCHLENGTH):
        """
        for each internal node, calculate the distance to each leaf,
        measured in branch length or internodes
        """
        if store is None:
            store = {}
        leaf2len = {}
        if self.children:
            for child in self.children:
                if measure == BRANCHLENGTH:
                    assert child.length is not None
                    dist = child.length
                elif measure == INTERNODES:
                    dist = 1
                else:
                    raise "InvalidMeasure"
                child.leaf_distances(store, measure)
                if child.istip:
                    leaf2len[child.label] = dist
                else:
                    for k, v in store[child].items():
                        leaf2len[k] = v + dist
        else:
            leaf2len[self] = {self.label: 0}
        store[self] = leaf2len
        return store

    def rootpath(self):
        n = self
        while 1:
            yield n
            if n.parent:
                n = n.parent
            else:
                break
            
    def subtree_mapping(self, labels, clean=False):
        """
        find the set of nodes in 'labels', and create a new tree
        representing the subtree connecting them.  nodes are assumed to be
        non-nested.

        return value is a mapping of old nodes to new nodes and vice versa.
        """
        d = {}
        oldtips = [ x for x in self.leaves() if x.label in labels ]
        for tip in oldtips:
            path = list(tip.rootpath())
            for node in path:
                if node not in d:
                    newnode = Node()
                    newnode.istip = node.istip
                    newnode.length = node.length
                    newnode.label = node.label
                    d[node] = newnode
                    d[newnode] = node
                else:
                    newnode = d[node]

                for child in node.children:
                    if child in d:
                        newchild = d[child]
                        if newchild not in newnode.children:
                            newnode.add_child(newchild)
        d["oldroot"] = self
        d["newroot"] = d[self]
        if clean:
            n = d["newroot"]
            while 1:
                if n.nchildren == 1:
                    oldnode = d[n]
                    del d[oldnode]; del d[n]
                    child = n.children[0]
                    child.parent = None
                    child.isroot = True
                    d["newroot"] = child
                    d["oldroot"] = d[child]
                    n = child
                else:
                    break
                    
            for tip in oldtips:
                newnode = d[tip]
                while 1:
                    newnode = newnode.parent
                    oldnode = d[newnode]
                    if newnode.nchildren == 1:
                        child = newnode.children[0]
                        if newnode.length:
                            child.length += newnode.length
                        newnode.remove_child(child)
                        if newnode.parent:
                            parent = newnode.parent
                            parent.remove_child(newnode)
                            parent.add_child(child)
                        del d[oldnode]; del d[newnode]
                    if not newnode.parent:
                        break
            
        return d

    def ultrametricize_dumbly(self):
        assert not self.istip
        d = self.leaf_distances()
        maxdist = max(d.values())
        for lf in self.leaves():
            x = d[lf.label]
            lf.length += maxdist - x

def node2size(node, d=None):
    "map node and descendants to number of descendant tips"
    if d is None:
        d = {}
    size = int(node.istip)
    if not node.istip:
        for child in node.children:
            node2size(child, d)
            size += d[child]
    d[node] = size
    return d

def reroot(oldroot, newroot):
    oldroot.isroot = False
    newroot.isroot = True
    v = []
    n = newroot
    while 1:
        v.append(n)
        if not n.parent: break
        n = n.parent
    #print [ x.label for x in v ]
    v.reverse()
    for i, cp in enumerate(v[:-1]):
        node = v[i+1]
        # node is current node; cp is current parent
        #print node.label, cp.label
        cp.remove_child(node)
        node.add_child(cp)
        cp.length = node.length
    return newroot

def bifid_tree(N):
    root = Node()
    root.isroot = True
    node = root
    for i in range(N):
        node.add_child(Node())
        node.add_child(Node())
        if i < (N-1):
            node.children[0].add_child(Node())
            node.children[0].add_child(Node())
        node = node.children[-1]
    i = 1
    for n in root.descendants():
        if not n.children:
            n.istip = True
            n.label = "SP%s" % i
        else:
            n.label = "IN%s" % i
        i += 1
    return root

def ultrametricize(node):
    import layout
    n2c = layout.calc_node_positions(node, 1.0, 1.0, scaled=False)
    for n in node.descendants():
        c = n2c[n]
        n.length = c.x - n2c[n.parent].x

if __name__ == "__main__":
    import newick, ascii, os
    ## from numpy import array
    ## #tree = newick.parse("(a,(b,(c,(d,e))));")
    ## f = os.path.expanduser("~/Projects/pedic-sympatry/matrices/")
    ## tree = eval(file(f+"garli-ml.tree").read())
    ## treespp = tree["species"]
    ## root = newick.parse(tree["newick"])
    ## spp = ['alaschanica', 'cheilanthifolia', 'dichotoma', 'kansuensis',
    ##        'oederi', 'plicata', 'przewalskii', 'remotiloba',
    ##        'rhinanthoides', 'roylei', 'rupicola', 'scolopax']
    ## print root.subtree_mapping(spp, clean=1)
    root = bifid_tree(5)
    ultrametricize(root)
    print newick.tostring(root)
    print ascii.render(root)
    #print node2tipsets(root)
    print root.mrca(("SP6","SP10", "SP18")).label
