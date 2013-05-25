#!/usr/bin/env python
import sys, math
from pprint import pprint
import scipy
from scipy import optimize
import phylo

LARGE = 10e10
PMAX = 100.0

def likelihood_de(params, model, tree):
    """
    demo function for optimizing dispersal and extinction rates, for
    use with scipy.optimize

    * returns negative log-likelihood of the tree, model, and data for
      dispersal and extinction rates given in params

    * assumes that the tree has only one model assigned to all branch
      segments
    """
    d, e = params
    for p in (d, e):
        if (p < 0) or (p > PMAX):
            return LARGE

    for p in params:
        if type(p) is complex:
            return LARGE
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    try:
        lh = tree.eval_likelihood()
        return -(scipy.log(lh))
    except:
        return LARGE

def likelihood_mp(params, model, tree):
    """
    """
    e = params[0]
    d = params[1]
        
    for p in params:
        if (p < 0) or (p > PMAX):
            return LARGE
        if type(p) is complex:
            return LARGE

    # The idea is that model.dp_array is the same shape as model.D,
    # and contains values that index free dispersal parameters, which
    # begin at index 2 in params (base extinction and dispersal are in
    # positions 0 and 1, respectively).
    if len(params) > 2:
        #v = scipy.array([1.0] + list(params[2:])) * d
        model.D = model.dp_array.choose(params[1:])
    else:
        model.setup_D(d)

    model.setup_E(e)
    model.setup_Q()
    try:
        lh = tree.eval_likelihood()
        #print params, -scipy.log(lh)
        return -(scipy.log(lh))
    except:
        return LARGE

def optimize_de(tree, model):
    """
    optimize dispersal and extinction rates
    model is an instance of RateModel
    """
    v = optimize.fmin_powell(
        likelihood_de, [0.01, 0.01], args=(model, tree),
        full_output=True, disp=1, callback=None
        )
    params, negloglikelihood = v[:2]
    if negloglikelihood == LARGE:
        raise Exception("ConvergenceError")
    dispersal, extinction = params
    return dispersal, extinction, negloglikelihood

def optimize_mp(tree, model):
    """
    optimize dispersal and extinction rates, multiparameter version
    model is an instance of DECModel
    """
    v = optimize.fmin_powell(
        likelihood_mp, model.params[:],
        args=(model, tree),
        full_output=True, disp=0
        )
    params, negloglikelihood = v[:2]
    if negloglikelihood == LARGE:
        raise Exception("ConvergenceError")

    return params, negloglikelihood

def ancsplit_likelihood_de(node, ancsplit, model, d, e):
    """
    calculate likelihood of ancsplit at node with dispersal and
    extinction rates d, e

    """
    #c1, c2 = node.children()
    c1, c2 = node.children
    seg1 = c1.segments[-1]; seg2 = c2.segments[-1]
    seg1.startdist = ancsplit.descdists[0]
    seg2.startdist = ancsplit.descdists[1]
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    lh = node.tree.eval_likelihood()
    seg1.startdist = None
    seg2.startdist = None
    return lh

def ancdist_likelihood_de(node, dist, model, d, e):
    """
    calculate likelihoods of ancsplits of dist at node with dispersal and
    extinction rates d, e (weighted average of split likelihoods)
    """
    #c1, c2 = node.children()
    c1, c2 = node.children
    seg1 = c1.segments[-1]; seg2 = c2.segments[-1]
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    v = []
    for ancsplit in model.iter_ancsplits(dist):
        seg1.startdist = ancsplit.descdists[0]
        seg2.startdist = ancsplit.descdists[1]
        lh = node.tree.eval_likelihood()
        #v.append(lh * ancsplit.weight)
        ancsplit.likelihood = lh
        v.append(ancsplit)
    seg1.startdist = None
    seg2.startdist = None
    return v

def ancdist_likelihood_mp(node, dist, model, params):
    """
    calculate likelihoods of ancsplits of dist at node with model params
    """
    #c1, c2 = node.children()
    c1, c2 = node.children
    seg1 = c1.segments[-1]; seg2 = c2.segments[-1]

    e = params[0]; d = params[1]
    if len(params) > 2:
        #v = scipy.array([1.0] + list(params[2:])) * d
        model.D = model.dp_array.choose(params[1:])
    else:
        model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    v = []
    for ancsplit in model.iter_ancsplits(dist):
        seg1.startdist = ancsplit.descdists[0]
        seg2.startdist = ancsplit.descdists[1]
        lh = node.tree.eval_likelihood()
        #v.append(lh * ancsplit.weight)
        ancsplit.likelihood = lh
        v.append(ancsplit)
    seg1.startdist = None
    seg2.startdist = None
    return v

def ancsplit_optimize_de(node, ancsplit, model):
    """
    optimize likelihood of ancsplit at node
    """
    #c1, c2 = node.children()
    c1, c2 = node.children
    seg1 = c1.segments[-1]; seg2 = c2.segments[-1]
    seg1.startdist = ancsplit.descdists[0]
    seg2.startdist = ancsplit.descdists[1]
    v = optimize.fmin_powell(
        likelihood_de, [0.01, 0.01], args=(model, node.tree),
        full_output=True,disp=0
        )
    seg1.startdist = None
    seg2.startdist = None
    params, negloglikelihood = v[:2]
    return params, negloglikelihood

def ancdist_optimize_de(node, ancdist, model):
    """
    optimize likelihood of ancdist at node
    """
    def f(params):
        d, e = params
        for p in (d, e):
            if (p < 0) or (p > PMAX):
                return LARGE
        try:
            return -scipy.log(ancdist_likelihood_de(node, ancdist, model, d, e))
        except:
            return LARGE
        
    v = optimize.fmin_powell(f, [0.01, 0.01], full_output=True,disp=0)
    params, negloglikelihood = v[:2]
    return params, negloglikelihood

def calculate_local_for_all_nodes(node, model):
    """
    optimize dispersal and extinction rates on tree with only one model

    * for each internal node, calculates likelihood and optimal
      dispersal and extinction rates for all split scenarios
    """
    tree = node.tree
    if not node.istip:
        #c1, c2 = node.children()
        c1, c2 = node.children
        calculate_local_for_all_nodes(c1, model)
        calculate_local_for_all_nodes(c2, model)
        print "Node %s --> nodes %s, %s" % (node.label, c1.label, c2.label)
        results = []
        tree.clear_startdist()
        for dist in model.dists:
            for ancsplit in model.iter_ancsplits(dist):
                d1, d2 = ancsplit.descdists
                c1.segments[-1].startdist = d1
                c2.segments[-1].startdist = d2
                v = optimize.fmin_powell(
                    likelihood_de, [0.1]*2, args=(model, tree),
                    full_output=True,disp=0
                    )
                params, opt = v[:2]
                #opt = scipy.log(exp(-opt) * weight) * -1
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                results.append((float(opt), tuple(params), ds1, ds2, root))
        results.sort()
        ptot = sum([ r[0] for r in results ])
        for opt, params, d1, d2, root in results:
            if opt > (results[0][0] + 2):
                break
            print "%s, %s; -lnL=%g (P=%g); d=%g; e=%g" % \
                  (d1, d2, opt, opt/ptotal, params[0], params[1])

def calculate_global_for_all_nodes(node, model, d, e, skip=[], _rv={}):
    if (not node.istip) and (not node.label in skip):
        #c1, c2 = node.children()
        c1, c2 = node.children
        calculate_global_for_all_nodes(c1, model, d, e)
        calculate_global_for_all_nodes(c2, model, d, e)
        print ", ".join([node.label, c1.label, c2.label])
        results = []
        #node.tree.clear_startdist()
        for disti, dist in model.enumerate_dists():
            S = set()
            for ancsplit in model.iter_ancsplits(dist):
                d1, d2 = ancsplit.descdists
                # calcuate likelihood of this split at this node
                lh = ancsplit_likelihood_de(node, ancsplit, model, d, e)
                
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(node.tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                try:
                    results.append((-(math.log(lh)), ds1, ds2, root))
                except:
                    pass

        results.sort()
        v = []
        for opt, ds1, ds2, root in results:
            if opt > (results[0][0] + 2):
                break
            v.append((opt, ds1, ds2))
            print "  -lnL %g, %s, %s" % (opt, ds1, ds2)
        _rv[node.label] = v

    return _rv


def ancsplits_mp(tree, model, params, nodelabels=None):
    node2splits = {}
    for node in [ x for x in tree.root.iternodes(phylo.PREORDER) \
                  if (not x.istip) ]:
        label = node.label
        if (nodelabels and (label in nodelabels)) or (not nodelabels):
            v = []
            for dist in [ d for d in model.dists[1:] \
                          if d not in node.excluded_dists ]:
                x = ancdist_likelihood_mp(node, dist, model, params)
                v.extend(x)
            v.sort(); v.reverse()
            node2splits[node] = v
    return node2splits
