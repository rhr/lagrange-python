import sys, os, nexus, scipy
from decmodel_mp import DECModel
from tree import Tree

VERSION = file(os.path.join(os.path.dirname(__file__),
                            "VERSION")).read().strip()

def eval_decmodel(s):
    d = eval(s)
    v = d["lagrange_version"] 
    if v != VERSION:
        print >> sys.stderr, "***Version mismatch: %s (expecting %s)***" \
              % (v, VERSION)
        print >> sys.stderr, "Things may not work as expected"
    for x in ("area_labels", "taxon_range_data", "newick_trees"):
        assert x in d, "required for analysis, but missing: %s" % x
    labels = d["area_labels"]
    nareas = len(labels)
    data = d["taxon_range_data"]
    maxareas = d.get("max_range_size") or len(labels)
    dists = d["ranges"]
    excluded = d["excluded_ranges"]
    dm = d["area_dispersal"]
    periods = d["dispersal_durations"]
    model = DECModel(nareas, labels, periods=periods, dists=dists)
    model.Dmask[:] = dm
    trees = d["newick_trees"]
    newicktree = trees[0]["newick"]
    nodelabels = trees[0]["included"]
    root_age = trees[0]["root_age"] or None
    tree = Tree(newicktree, periods=periods, root_age=root_age)
    tree.set_default_model(model)
    tree.set_tip_conditionals(data)
    base_rates = d["base_rates"]
    if base_rates == "__estimate__":
        pass
    else:
        base_rates = (d["base_rates"]["dispersal"], d["base_rates"]["extinction"])
    return model, tree, data, nodelabels, base_rates
    
def parse_upper_triangle(s):
    """
    assumes first row is labels,
    and each subsequent row's first token is a label
    """
    try:
        import graph
    except:
        print >> sys.stderr, "igraph library not available; see igraph.sf.net"
        return
    lines = s.split("\n")
    labels = lines.pop(0).split()
    nareas = len(labels)
    dm = scipy.ones((nareas, nareas))
    edges = []
    for line in lines:
        row = line.split()
        t = row.pop(0)
        i = labels.index(t)
        for j, v in enumerate(row):
            v = v.strip()
            if (i != j) and v:
                dm[i,j] = float(v); dm[j,i] = float(v)
                edges.append((i,j))
    return dm, graph.AreaGraph(edges=edges)

def parse_lower_triangle(s, labels=None):
    try:
        import graph
    except:
        print >> sys.stderr, "igraph library not available; see igraph.sf.net"
        return
    tokens = [ x.strip().split() for x in s.strip().split("\n") ]
    firstrow = tokens[0]
    if (len(firstrow) == 1 and firstrow[0]=='-') or \
           (labels and firstrow[0] in labels \
            and len(firstrow) == 2 and firstrow[1] == '-'):
        tokens = tokens[1:]
    nareas = len(tokens)+1
    #print "nareas", nareas; import sys; sys.exit()
    dm = scipy.ones((nareas, nareas))
    
    edges = []
    for ai, v in enumerate(tokens):
        if v[-1] == "-":
            v = v[:-1]
            ai += 1
        if labels:
            label = v.pop(0)
            ai = labels.index(label)
        for i, x in enumerate(v):
            x = float(x)
            dm[ai,i] = x; dm[i,ai] = x
            if x > 0:
                edges.append((ai,i))
    return dm, graph.AreaGraph(edges=edges)

def parse_nexus(infile):
    return nexus.Nexus(infile)

def parse_aln(infile):
    ntax = 0; nareas = 0
    taxon2dist = {}
    for line in [ x.strip() for x in infile if x.strip() ]:
        if line.startswith("#"):
            continue
        try:
            ntax, nareas = map(int, line.split())
            continue
        except:
            pass
        taxon, s = line.split()
        dist = tuple([ bool(int(x)) for x in s ])
        taxon2dist[taxon] = dist
    return taxon2dist

def parse_labeldata(infile, labels, labelsep=None):
    if type(labels) is str:
        labels = list(labels)
    taxon2dist = {}; nareas = len(labels)
    for line in [ x.strip() for x in infile if x.strip() ]:
        #print line
        if line.startswith("#"):
            continue
        taxon, areas = line.split()
        if labelsep:
            areas = areas.split(labelsep)
        dist = [False]*nareas
        for label in areas:
            dist[labels.index(label)] = True
        assert taxon not in taxon2dist, "duplicate taxon %s" % taxon
        taxon2dist[taxon] = tuple(dist)
    return taxon2dist

def parse_matrix(infile):
    labels = []
    taxon2dist = {}
    for line in [ x.strip() for x in infile if x.strip() ]:
        if line.startswith("#"):
            continue
        tokens = line.split()

        if not labels:
            if len(tokens) == 1:
                labels = list(tokens[0])
            else:
                labels = tokens
            nareas = len(labels)
            continue

        else:
            if len(tokens) == 2:
                tokens = [tokens[0]] + list(tokens[1])

        taxon = tokens.pop(0)
        assert taxon not in taxon2dist, \
               "Duplicate taxon '%s' in data file" % taxon
        dist = [False]*nareas
        assert len(tokens) == nareas
        for i, t in enumerate(tokens):
            dist[i] = bool(int(t))
        assert taxon not in taxon2dist, "duplicate taxon %s" % taxon
        taxon2dist[taxon] = tuple(dist)
    return labels, taxon2dist

def parse_matrix2(infile):
    """
    return area_labels and data (dictionary mapping taxa to their
    observed ranges)
    """
    labels = []
    taxon2dist = {}
    for line in [ x.strip() for x in infile if x.strip() ]:
        if line.startswith("#"):
            continue
        tokens = line.split()

        if not labels:
            if len(tokens) == 1:
                labels = list(tokens[0])
            else:
                labels = tokens
            nareas = len(labels)
            continue

        else:
            if len(tokens) == 2:
                tokens = [tokens[0]] + list(tokens[1])

        taxon = tokens.pop(0)
        assert taxon not in taxon2dist, \
               "Duplicate taxon '%s' in data file" % taxon
        dist = []
        assert len(tokens) == nareas
        for i, t in enumerate(tokens):
            if bool(int(t)):
                dist.append(i)
        assert taxon not in taxon2dist, "duplicate taxon %s" % taxon
        taxon2dist[taxon] = tuple(dist)
    return labels, taxon2dist

def parse_matrix3(infile):
    "return area_labels, taxa, and data"
    labels = []
    taxa = []
    taxon2dist = {}
    for line in [ x.strip() for x in infile if x.strip() ]:
        if line.startswith("#"):
            continue
        tokens = line.split()

        if not labels:
            if len(tokens) == 1:
                labels = list(tokens[0])
            else:
                labels = tokens
            nareas = len(labels)
            continue

        else:
            if len(tokens) == 2:
                tokens = [tokens[0]] + list(tokens[1])

        taxon = tokens.pop(0)
        assert taxon not in taxon2dist, \
               "Duplicate taxon '%s' in data file" % taxon
        taxa.append(taxon)
        dist = []
        assert len(tokens) == nareas
        for i, t in enumerate(tokens):
            if bool(int(t)):
                dist.append(i)
        assert taxon not in taxon2dist, "duplicate taxon %s" % taxon
        taxon2dist[taxon] = tuple(dist)
    return labels, taxa, taxon2dist
        

if __name__ == "__main__":
    #print parse_aln(file("../psychotria.aln"))
    #print file("../psychotria.data").read()
    #print parse_labeldata(file("../psychotria.data"), "KOMH")
    print parse_matrix(file("../psychotria.matrix"))
