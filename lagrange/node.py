"""
    new node class
"""
PREORDER = 0
POSTORDER = 1

class Node:
    def __init__(self, parent = None, childs = None, isroot = 0,
                 label = "", length = None, istip = 1):
        self.parent = parent
        self.childs = childs
        self.isroot = isroot
        self.label = label
        self.length = length
        self.istip = istip
        self.distanceToTip = 0

    def unlink(self):
        """Remove references to next, back, and data, to let refcounts
        go to zero"""
        del self.childs; del self.parent

    def fnodes(self):
        """Returns a list of the Fnodes linked to self, including self"""
        nodes = []
        n = self
        if self.parent != None:
            nodes.append(self.parent)
        if self.childs != None:
            for i in range(len(self.childs)):
                nodes.append(self.childs[i])
        return nodes

    def iterchildren(self):
        """Returns the immediate descendants of this node"""
        if not self.istip:
            for i in range(len(childs)):
                yield childs[i]
                
    def children(self):
        """Returns the immediate descendants of this node"""
        if self.childs == None:
            return []
        return self.childs

    def iterdescendants(self, order=PREORDER):
        """Returns a list of all descendants of this node."""
        if order == PREORDER:
            yield self
        if self.children() != None:
            for child in self.children():
                for d in child.descendants(order):
                    yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER):
        """Returns a list of all descendants of this node."""
        return list(self.iterdescendants(order))
                
    def iterleaves(self):
        """Returns a list of leaf nodes that are descendant from this
        node."""
        for child in self.children():
            if child.istip:
                yield child
            else:
                for leaf in child.leaves():
                    yield leaf

    def leaves(self):
        """Returns a list of leaf nodes that are descendant from this
        node."""
        return list(self.iterleaves())

    def add_child(self, node):
        """Adds node to self's children"""
        self.istip = 0
        if self.childs == None:
            self.childs = []
        if self.childs.count(node) < 1:
            self.childs.append(node)
            node.parent = self
            return True
        else:
            return False

    def remove_child(self, node):
        """Adds node to self's children"""
        if self.childs.count(node) > 0:
            self.childs.remove(node)
            node.parent = None
            return true
        else:
            return false

    def prune(self):
        # use to return last and next, don't know what to do now
        # don't do with root
        # assumes bisecting tree, need to fix
        assert self.parent
        parent = self.parent
        if parent != self.root:
            child = None
            for i in range(len(parent(children()))):
                if parent.children[i] != node:
                    child = parent.children[i]
            pparent = parent.parent
            parent.remove_child(node)
            parent.remove_child(child)
            pparent.remove_child(parent)
            pparent.add_child(child)
            child.length = child.length+parent.length
        else:
            child = None
            for i in range(len(parent(children()))):
                if parent.children[i] != node:
                    child = parent.children[i]
            parent.remove_child(node)
            child.parent = None

    def bisect(self):
        """Bisect the branch between self and self.back with a new
        internal node"""
        assert self.parent
        parent = self.parent
        n = InternalNode()
        # join the new internal node to self and self.back
        n.parent = parent; self.parent = n
        for i in range(len(parent.children())):
            if parent.children()[i] == self:
                parent.children()[i] = n
        length = self.length
        if length:
            half = length*0.5
            self.length = half; n.length = half
        return n

    def has_descendant(self, node):
        """check of node is descendant of self"""
        if self.istip: return 0
        flag = 0
        if self.childs.count(node) > 0:
            flag = 1
        return flag

    def print_leaves(self):
        for leaf in self.leaves():
            print leaf.label,
        print

def polarize(node, parent=None):
    node.parent = parent
    if not node.istip:
        for child in node.children():
            polarize(child, node)

def InternalNode(isroot=0):
    """return an internal node"""
    node = Node(isroot=isroot)
    return node        

def getDistanceFromTip(node):
    cur = node
    curh = 0.0
    while cur != None:
        curh += cur.length
        cur = cur.parent
    return curh
