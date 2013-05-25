import scipy
import phylo, newick, ascii, optimize

def log(s, outfile, tee=False):
    print >> outfile, s
    if tee:
        print s

def summarize_splits(splits, weighted=True):
    rows = []
    if weighted:
        v = [ (x.likelihood * x.weight, x) for x in splits ]
    else:
        v = [ (x.likelihood, x) for x in splits ]
    ptot = sum([ x[0] for x in v ])
    v.sort(); v.reverse()
    opt = scipy.log(v[0][0])

    rows.append(["split", "lnL", "Rel.Prob"])
    sumprob = 0.0
    for L, split in v:
        lnL = scipy.log(L)
        relprob = (L/ptot)
        if sumprob < 0.95:
        #if (opt - lnL) < 2:
            rows.append([str(split), "%.4g" % lnL, "%.4g" % relprob])
        sumprob += relprob
    widths = []
    for i in range(3):
        w = max([ len(x[i]) for x in rows ])
        for x in rows:
            x[i] = x[i].ljust(w)
    return [ "  ".join(x) for x in rows ]

def ascii_tree(outfile, tree, model, data=None, scaled=False, minwidth=80,
               tee=False):
    d = None
    if data:
        d = dict([ (x[0], model.dist2label(x[1])) for x in data.items() \
                   if x[1] in model.dists ])
    s = "\n".join(
        ["Newick tree with interior nodes labeled:",
         newick.to_string(tree.root)+";", "\n",
         "Cladogram (branch lengths not to scale):",
         ascii.render(tree.root, scaled=scaled, minwidth=80, data=d),
         "\n\n"]
        )
    log(s, outfile, tee)

def ancsplits(outfile, tree, model, dispersal, extinction,
              nodelabels=None, verbose=True, tee=False):
    # nodelabels is assumed to be a list of strings
    if verbose:
        s = """
Ancestral range subdivision/inheritance scenarios ('splits') at
internal nodes.

* Split format: [left|right], where 'left' and 'right' are the ranges
  inherited by each descendant branch (on the printed tree, 'left' is
  the upper branch, and 'right' the lower branch).

* Only splits within 2 log-likelihood units of the maximum for each
  node are shown.  'Rel.Prob' is the relative probability (fraction of
  the global likelihood) of a split.
"""
        log(s, outfile, tee)
    for node in [ x for x in tree.root.iternodes(phylo.PREORDER) \
                  if (not x.istip) ]:
#              if (not x.istip) and (x.parent) ]:
        #c1, c2 = node.children
        label = node.label
        if (nodelabels and (label in nodelabels)) or (not nodelabels):
            #print >> outfile, "At node %s:" % label
            log("At node %s:" % label, outfile, tee)
            v = []
            for dist in [ d for d in model.dists[1:] \
                          if d not in node.excluded_dists ]:
                x = optimize.ancdist_likelihood_de(
                    node, dist, model, dispersal, extinction
                    )
                v.extend(x)
            for line in summarize_splits(v, weighted=False):
                #print >> outfile, "  ", line
                log("   %s" % line, outfile, tee)
            #print >> outfile, ""
            log("", outfile, tee)
            outfile.flush()

def ancsplits_mp(outfile, tree, model, params,
                 nodelabels=None, verbose=True, tee=False):
    # nodelabels is assumed to be a list of strings
    if verbose:
        s = """
Ancestral range subdivision/inheritance scenarios ('splits') at
internal nodes.

* Split format: [left|right], where 'left' and 'right' are the ranges
  inherited by each descendant branch (on the printed tree, 'left' is
  the upper branch, and 'right' the lower branch).

* Only splits within 2 log-likelihood units of the maximum for each
  node are shown.  'Rel.Prob' is the relative probability (fraction of
  the global likelihood) of a split.
"""
        log(s, outfile, tee)
    for node in [ x for x in tree.root.iternodes(phylo.PREORDER) \
                  if (not x.istip) ]:
#              if (not x.istip) and (x.parent) ]:
        #c1, c2 = node.children
        label = node.label
        if (nodelabels and (label in nodelabels)) or (not nodelabels):
            #print >> outfile, "At node %s:" % label
            log("At node %s:" % label, outfile, tee)
            v = []
            for dist in [ d for d in model.dists[1:] \
                          if d not in node.excluded_dists ]:
                x = optimize.ancdist_likelihood_mp(
                    node, dist, model, params
                    )
                v.extend(x)
            for line in summarize_splits(v, weighted=False):
                log("   %s" % line, outfile, tee)
            log("", outfile, tee)
            outfile.flush()

def optimize_dispersal_extinction(outfile, tree, model, tee=False):
    d, e, nlnL = optimize.optimize_de(tree, model)
    #print >> outfile, "Global ML at root node:\n  -lnL = %.4g" % nlnL
    log("Global ML at root node:\n  -lnL = %.4g" % nlnL, outfile, tee)
    #print >> outfile, "  dispersal = %.4g\n  extinction = %.4g" % (d, e)
    log("  dispersal = %.4g\n  extinction = %.4g" % (d, e), outfile, tee)
    return d, e

def optimize_dispersal_extinction_mp(outfile, tree, model, tee=False):
    params, nlnL = optimize.optimize_mp(tree, model)
    e, d = params
    #print >> outfile, "Global ML at root node:\n  -lnL = %.4g" % nlnL
    log("Global ML at root node:\n  -lnL = %.4g" % nlnL, outfile, tee)
    #print >> outfile, "  dispersal = %.4g\n  extinction = %.4g" % (d, e)
    log("  dispersal = %.4g\n  extinction = %.4g" % (d, e), outfile, tee)
    return d, e
