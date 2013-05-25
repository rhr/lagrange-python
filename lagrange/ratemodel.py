#!/usr/bin/env python
"""
conventions:

dist = shorthand for distribution = binary vector indicating
       presence/absence in a vector of areas

"""
import sys, sets, string, time
import scipy
import rates
import nchoosem

class RateModel:
    """
    Model of dispersal and local extinction between discrete
    geographic areas.
    """
    def __init__(self,
                 nareas,         # number of areas
                 labels=None,    # ordered sequence of area labels
                 periods=None,   # list of durations, ordered from
                                 # present to past, corresponding to
                                 # time-slices of dispersal and
                                 # extinction parameters; sum(periods)
                                 # must exceed the root age of any
                                 # tree that uses the model
                 dists=None,
                 Dmask=None
                 ):
        """
        initialize variables
        """
        self.nareas = nareas
        self.arange = range(nareas)
        if not labels:
            if nareas < 26:
                self.labels = string.uppercase[:nareas]
            else:
                self.labels = [ "A%d" % i for i in self.arange ]
        else:
            assert len(labels) == nareas, \
                   "Mismatch between number of areas and number of area labels"
            self.labels = labels
        if type(self.labels) == str or \
               max([ len(x) for x in labels ]) == 1:
            self.labelsep = ""
        else:
            self.labelsep = "+"
        self.label2i = dict([ (label, i) for i, label \
                              in enumerate(self.labels) ])
        self.periods = periods or []
        self.nperiods = len(self.periods) or 1
        self.prange = range(self.nperiods)
        # set up self.dists, the list of dists corresponding to the
        # number of areas; does not include the empty dist (vector of
        # all zeros) (cf. RateModelGE)
        self.setup_dists(dists)

        pi = 1.0/self.ndists
        self.dist_priors = [ pi for x in self.distrange ]

        #print self.dist_priors

        # instantiate a matrix for each period
        # of nareas x nareas
        self.Dmask = scipy.ones((self.nperiods, nareas, nareas))
        if Dmask is not None:
            # make sure default Dmask has correct dimensions
            assert Dmask.shape == (nareas, nareas)
            self.Dmask[:] = Dmask
        default_d = 0.01
        default_e = 0.01
        self.setup_D(default_d)
        self.setup_E(default_e)
        #t = time.time()
        self.setup_Q()
        #print "%s for setup_Q" % (time.time() - t)

    def label2dist(self, label):
        assert label in self.diststrings, \
               "Range '%s' not defined in model" % label
        v = self.dists[self.diststrings.index(label)]
        return v


    def dist2label(self, dist):
        return self.diststrings[self.dist2i[dist]]

    def setup_dists(self, dists=None):
        # results in somthing like 
        # [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        # for three areas
        if dists is None:
            self.dists = [ dist for dist in \
                           nchoosem.iterate_all_bv(self.nareas) ]
        else:
            self.dists = dists[:]
        self.diststrings = []
        for d in self.dists:
            self.diststrings.append(self.labelsep.join(
                [ self.labels[i] for i in self.arange if d[i] ]
                ))
                
        self.ndists = len(self.dists)
        self.dist2i = dict([
            (dist, i) for i, dist in enumerate(self.dists)
            ])
        self.distrange = range(self.ndists)

    def setup_D(self, d):
        nareas = self.nareas
        self.D = scipy.ones((self.nperiods, nareas, nareas)) * d
        #for period, duration in enumerate(self.periods):
        for period in self.prange:
            for i in self.arange:
                self.D[period,i,i] = 0.0
        self.D *= self.Dmask

#these are both wrong - should take area index as param, not dists
##     def set_D_cell(self, period, fromdist, todist, d):
##         i = self.dist2i[fromdist]
##         j = self.dist2i[todist]
##         self.D[period, i, j] = d

##     def set_Dmask_cell(self, period, fromdist, todist, d, sym=False):
##         i = self.dist2i[fromdist]
##         j = self.dist2i[todist]
##         self.Dmask[period, i, j] = d
##         if sym:
##             self.Dmask[period, j, i] = d

    def set_Dmask_cell(self, period, a1, a2, d, sym=False):
        if a1 in self.labels:
            i = self.label2i[a1]
        else:
            i = a1
        if a2 in self.labels:
            j = self.label2i[a2]
        else:
            j = a2
        self.Dmask[period, i, j] = d
        if sym:
            self.Dmask[period, j, i] = d

    def setup_E(self, e):
        self.E = scipy.ones((self.nperiods, self.nareas)) * e

    def setup_Q(self):
        t = 0.0
        self.Q = scipy.zeros((self.nperiods, self.ndists, self.ndists))
        #D = self.D * self.Dmask
        logical_xor = scipy.logical_xor
        nonzero = scipy.nonzero
        for p in range(self.nperiods):
            for i, d1 in enumerate(self.dists):
                s1 = sum(d1)
                if s1 > 0:
                    for j, d2 in enumerate(self.dists):
                        s2 = sum(d2)
                        xor = logical_xor(d1, d2)
                        # only consider transitions between dists that are
                        # 1 step (dispersal or extinction) apart
                        if sum(xor) == 1:
                            dest = nonzero(xor)[0]
                            #prior = self.dist_priors[i]
                            if s1 < s2: # dispersal
                                rate = 0.0
                                #t1 = time.time()
                                for src in nonzero(d1)[0]:
                                    rate += self.D[p,src,dest] * \
                                            self.Dmask[p,src,dest]
                                #t += time.time() - t1
                                # for each area in d1, add rate of
                                # dispersal to dest
                            else: # extinction
                                rate = self.E[p,dest]
                            #self.Q[i][j] = (prior * rate)
                            try:
                                self.Q[p,i,j] = rate
                            except ValueError:
                                self.Q[p,i,j] = abs(rate)
            self.set_Qdiag(p)
        #print "t = %g" % t

    def set_Qdiag(self, period):
        for i in self.distrange:
            self.Q[period,i,i] = (sum(self.Q[period,i,:]) - \
                                  self.Q[period,i,i]) * -1.0

    def P(self, period, t):
        """
        return P, the matrix of dist-to-dist transition probabilities,
        from the model's rate matrix (Q) over a time duration (t)
        """
        p = rates.Q2P(self.Q[period], t)
        # filter out impossible dists
        for i, dist in self.enumerate_dists():
            for area, present in enumerate(dist):
                if present and (sum(self.Dmask[period,area,:])+\
                                sum(self.Dmask[period,:,area]) == 0):
                    p[i,:] *= 0.0
                    break
        return p

    def Pold(self, period, t):
        """
        return P, the matrix of dist-to-dist transition probabilities,
        from the model's rate matrix (Q) over a time duration (t)
        """
        p = rates.Q2P_old(self.Q[period], t)
        # filter out impossible dists
        for i, dist in self.enumerate_dists():
            for area, present in enumerate(dist):
                if present and (sum(self.Dmask[period,area,:])+\
                                sum(self.Dmask[period,:,area]) == 0):
                    p[i,:] *= 0.0
                    break
        return p

    def Q_repr(self, period):
        lines = []
        widths = [ max([ max((len("%g" % x), self.nareas)) \
                         for x in self.Q[period,:,col] ]) \
                   for col in self.distrange ]
        lines.append("  ".join(
            [ ("".join([ str(a) for a in d ])).rjust(widths[col]) \
              for col, d in enumerate(self.dists) ]))
        lines.append("".join(["-"] * len(lines[0])))
                     
        for row in self.distrange:
            lines.append(
                "  ".join([ ("%g" % x).rjust(widths[col]) \
                            for col, x in enumerate(self.Q[period,row,:]) ])
                )

        w = max(len(r"From\To"), self.nareas)
        for i, line in enumerate(lines):
            if i == 0:
                lines[i] = "%s %s" % (r"From\To".rjust(w), line)
            elif i == 1:
                lines[i] = "%s%s" % ("".join(["-"]*(w+1)), line)
            else:
                d = self.dists[i-2]
                lines[i] = "%s|%s" % ("".join([ str(a) for a in d ]).ljust(w),
                                      line)

        return "\n".join(lines)

    def P_repr(self, period, t):
        P = self.P(period, t)#rates.Q2P(self.Q[period], t)
        lines = []
        widths = [ max([ max((len("%g" % x), self.nareas)) \
                         for x in P[:,col] ]) \
                   for col in self.distrange ]
        lines.append("  ".join(
            [ ("".join([ str(a) for a in d ])).rjust(widths[col]) \
              for col, d in enumerate(self.dists) ]))
        lines.append("".join(["-"] * len(lines[0])))
                     
        for row in self.distrange:
            lines.append(
                "  ".join([ ("%g" % x).rjust(widths[col]) \
                            for col, x in enumerate(P[row,:]) ])
                )

        w = max(len(r"From\To"), self.nareas)
        for i, line in enumerate(lines):
            if i == 0:
                lines[i] = "%s %s" % (r"From\To".rjust(w), line)
            elif i == 1:
                lines[i] = "%s%s" % ("".join(["-"]*(w+1)), line)
            else:
                d = self.dists[i-2]
                lines[i] = "%s|%s" % ("".join([ str(a) for a in d ]).ljust(w),
                                      line)

        return "\n".join(lines)

    def enumerate_dists(self):
        "enumerate non-empty dists"
        return [ (i, d) for i, d in zip(self.distrange, self.dists) \
                 if sum(d) ]

    def iter_dist_splits(self, dist):
        assert dist in self.dists
        if sum(dist) == 1:
            yield (dist, dist)
        else:
            for i in scipy.nonzero(dist)[0]:
                x = scipy.zeros((len(dist),), dtype="i")
                x[i] = 1
                x = tuple(x)
                if x in self.dists:
                    yield (x, dist)
                    yield (dist, x)
                    y = tuple(scipy.array(scipy.logical_xor(dist, x),
                                          dtype="i"))
                    if y in self.dists:
                        yield (x, y)
                        if sum(y) > 1:
                            yield (y, x)
                    
    def iter_ancsplits(self, dist):
        splits = [ s for s in self.iter_dist_splits(dist) ]
        if splits:
            nsplits = len(splits)
            weight = (1.0/nsplits)
            for split in splits:
                as = Ancsplit(self, dist, split, weight=weight)
                yield as

    def remove_dist(self, dist):
        t = type(dist)
        if t is str:
            self.dists.remove(self.label2dist(dist))
        else:
            self.dists.remove(dist)
        self.setup_dists(self.dists)

class RateModelGE(RateModel):
    """
    Subclass of RateModel that allows global extinction
    """
        
    def setup_dists(self, dists=None):
        if dists is None:
            self.dists = [ dist for dist in \
                           nchoosem.iterate_all_bv2(self.nareas) ]
        else:
            self.dists = dists[:]
            emptydist = tuple([0]*self.nareas)
            if emptydist not in dists:
                self.dists.insert(0, emptydist)

        self.diststrings = []
        for d in self.dists:
            self.diststrings.append(self.labelsep.join(
                [ self.labels[i] for i in self.arange if d[i] ]
                ))
                
        self.ndists = len(self.dists)
        self.dist2i = dict([
            (dist, i) for i, dist in enumerate(self.dists)
            ])
        self.distrange = range(self.ndists)

## print RateModelGE(3).enumerate_dists()
## sys.exit()

class RateModelGE2(RateModelGE):
    def setup_dists(self, dists=None):
        if dists is None:
            self.dists = [ dist for dist in \
                           nchoosem.iterate_all_idx(self.nareas) ]
        else:
            self.dists = dists[:]

        emptydist = tuple()
        if emptydist not in self.dists:
            self.dists.insert(0, emptydist)

        self.diststrings = []
        for d in self.dists:
            self.diststrings.append(self.labelsep.join(
                [ self.labels[i] for i in d ]
                ))
                
        self.ndists = len(self.dists)
        self.dist2i = dict([
            (dist, i) for i, dist in enumerate(self.dists)
            ])
        self.distrange = range(self.ndists)

        self.dist_dispersals = set()
        self.dist_extinctions = set()
        for i, d1 in enumerate(self.dists):
            s1 = len(d1)
            if s1 > 0:
                for j, d2 in enumerate(self.dists):
                    s2 = len(d2)
                    if s1 == s2+1:
                        diff = set(d1) - set(d2)
                        if len(diff) == 1:
                            dest = tuple(diff)[0]
                            self.dist_extinctions.add((i,j,dest))
                            self.dist_dispersals.add((j,i,dest))
                    elif s1 == s2-1:
                        diff = set(d2) - set(d1)
                        if len(diff) == 1:
                            dest = tuple(diff)[0]
                            self.dist_dispersals.add((i,j,dest))
                            self.dist_extinctions.add((j,i,dest))


    def setup_Q(self):
        self.Pdict = {}
        self.Q = scipy.zeros((self.nperiods, self.ndists, self.ndists))
        #D = self.D * self.Dmask
        for p in range(self.nperiods):
            for i, j, dest in self.dist_dispersals:
                d1 = self.dists[i]
                rate = 0.0
                for src in d1:
                    rate += self.D[p,src,dest] * self.Dmask[p,src,dest]
                try:
                    self.Q[p,i,j] = rate
                except ValueError:
                    self.Q[p,i,j] = abs(rate)

            for i, j, dest in self.dist_extinctions:
                rate = self.E[p,dest]
                try:
                    self.Q[p,i,j] = rate
                except ValueError:
                    self.Q[p,i,j] = abs(rate)
                
            self.set_Qdiag(p)


    def enumerate_dists(self):
        "enumerate non-empty dists"
        return [ (i, d) for i, d in zip(self.distrange, self.dists) \
                 if len(d) ]

    def P(self, period, t):
        """
        return P, the matrix of dist-to-dist transition probabilities,
        from the model's rate matrix (Q) over a time duration (t)
        """
        k = (period, t) 
        if k in self.Pdict:
            return self.Pdict[k]
        p = rates.Q2P(self.Q[period], t)
        # filter out impossible dists
        for i, dist in self.enumerate_dists():
            for area in dist:
                if (sum(self.Dmask[period,area,:])+\
                    sum(self.Dmask[period,:,area])) == 0:
                    p[i,:] *= 0.0
                    break
        self.Pdict[k] = p
        return p

    def iter_dist_splits(self, dist):
        assert dist in self.dists
        if len(dist) == 1:
            yield (dist, dist)
        else:
            for i in dist:
                x = (i,)
                if x in self.dists:
                    yield (x, dist)
                    yield (dist, x)
                    ldist = list(dist)
                    ldist.remove(i)
                    ldist = tuple(ldist)
                    if ldist in self.dists:
                        yield (x, ldist)
                        if len(ldist) > 1:
                            yield (ldist, x)
                    
    def iter_ancsplits(self, dist):
        splits = [ s for s in self.iter_dist_splits(dist) ]
        if splits:
            nsplits = len(splits)
            weight = (1.0/nsplits)
            for split in splits:
                as = Ancsplit(self, dist, split, weight=weight)
                yield as


class Ancsplit:
    """
    convenience class for encapsulating an ancestor range splitting
    into descendant ranges
    """
    def __init__(self, model, ancdist, descdists, weight=None, likelihood=None):
        self.model = model
        self.ancdist = ancdist
        self.descdists = descdists
        self.weight = weight
        self.likelihood = likelihood

    def __repr__(self):
        d1, d2 = map(self.model.dist2label, self.descdists)
        lh = self.likelihood
        if lh: lh = "%.3g" % lh
        w = self.weight
        if w: w = "%.3g" % w
        return "[%s|%s]" % (d1, d2)

#
# utility functions
#

def iter_dist_splits(dist):
    if sum(dist) == 1:
        yield (dist, dist)
    else:
        for i in scipy.nonzero(dist)[0]:
            x = scipy.zeros((len(dist),), dtype="i")
            x[i] = 1
            x = tuple(x)
            yield (x, dist)
            yield (dist, x)
            y = tuple(scipy.array(scipy.logical_xor(dist, x), dtype="i"))
            yield (x, y)
            if sum(y) > 1:
                yield (y, x)

## d = (1,1,0,0)
## for x in iter_dist_splits(d):
##     print x
## sys.exit()

def dist_splits(dist):
    return sets.Set([ s for s in iter_dist_splits(dist) ])

## def iter_ancsplits(dist):
##     splits = [ s for s in iter_dist_splits(dist) ]
##     nsplits = len(splits)
##     weight = (1.0/nsplits)
##     for split in splits:
##         as = Ancsplit(dist, split, weight=weight)
##         yield as

## d = (1,1,0,0)
## for x in iter_ancsplits(d):
##     print x
## sys.exit()

def iter_dist_splits_weighted(dist):
    s = sum(dist)
    if s == 1:
        yield (dist, dist), 1.0
    else:
        wt = 1.0/(s*4)
        for sp in iter_dist_splits(dist):
            yield sp, wt

## d = (1,1,0)
## for x in iter_dist_splits_weighted(d):
##     print x
## print
## sys.exit()

def test_conditionals(distconds, seglens, model):
    """
    small test function to make sure evaluating likelihoods along
    segmented branches works properly
    """
    distrange = model.distrange
    for p, seglen in enumerate(seglens):
        print "distconds", distconds
        P = model.P(p, seglen)
        v = scipy.zeros((model.ndists,))
        for i in distrange:
            # P[i] is the vector of probabilities of going from dist i
            # to all other dists
            v[i] = sum(distconds * P[i])
        distconds = v
    return distconds

## def conditionals(node):
##     """
##     calculate the conditional likelihoods of ranges at the start of a
##     branch given the likelihoods of states at its end
##     """

##     # conditional likelihoods of dists at the end of the branch: for
##     # tips, this is a zeroed vector with 1.0 at the observed range
##     distconds = node.segments[0].dist_conditionals
##     #print node.label, node.segments, node.length
##     #for model, time, startdist in node.segments:
##     for seg in node.segments:
##         seg.dist_conditionals = distconds
##         model = seg.model
##         P = model.P(seg.period, seg.duration)
##         v = scipy.zeros((model.ndists,))
##         if not seg.distrange:
##             if seg.startdist:
##                 seg.distrange = (model.dist2i[seg.startdist],)
##             elif seg.fossils:
##                 # seg.fossils is a list of area indices
##                 print "fossils at node %s: %s" % (node.label, seg.fossils)
##                 num_fossil_areas = len(seg.fossils)
##                 distrange = list(model.distrange)
##                 excluded = []
##                 for i, dist in model.enumerate_dists():
##                     if sum(scipy.take(dist, seg.fossils)) != num_fossil_areas:
##                         distrange.remove(i)
##                         excluded.append(i)
##                 seg.distrange = distrange
##                 print "ranges excluded: %s" % excluded
##                 sys.exit()
##             else:
##                 seg.distrange = model.distrange

##         for i in seg.distrange:
##             # P[i] is the vector of probabilities of going from dist i
##             # to all other dists
##             v[i] = sum(distconds * P[i])
##         distconds = v
##     return distconds

## def test_conditionals(model, period, duration, distconds):
##     P = model.P(period, duration)
##     v = scipy.zeros((model.ndists,))
##     distrange = model.distrange
##     for i in distrange:
##         # P[i] is the vector of probabilities of going from dist i
##         # to all other dists
##         v[i] = sum(distconds * P[i])
##     return v

## def ancdist_conditional_lh(node):
##     """
##     recursive calculation of fractional likelihoods for dists at
##     internal nodes in a tree
##     """
##     if not node.istip:
##         c1, c2 = node.children()
##         if node.parent:
##             model = node.segments[0].model
##         else:
##             model = node.model
##         dist2i = model.dist2i

##         ancdist_conditional_lh2(c1)
##         ancdist_conditional_lh2(c2)

##         v1 = conditionals(c1)
##         v2 = conditionals(c2)

##         distconds = scipy.zeros((model.ndists,))
##         ancsplits = []
##         for distidx, dist in model.enumerate_dists():
##             lh = 0.0
##             pi = model.dist_priors[distidx]
## #            if dist not in node.excluded_dists:
##             for ancsplit in model.iter_ancsplits(dist):
##                 d1, d2 = ancsplit.descdists
##                 lh_part = (v1[dist2i[d1]] * v2[dist2i[d2]])
##                 lh += (lh_part * ancsplit.weight)
##                 ancsplit.likelihood = lh_part
##                 ancsplits.append(ancsplit)
##             #distconds[distidx] = lh
##             distconds[distidx] = lh * pi
##         node.ancsplits = ancsplits

##     else:
##         distconds = node.segments[0].dist_conditionals
    
##     if node.parent:
##         node.segments[0].dist_conditionals = distconds
##     else:
##         node.dist_conditionals = distconds

def nondiag_indices(m):
    """
    m is a square array - iterate over the (i,j) pairs indexing the
    non-diagonal cells
    """
    ind = 0
    N = len(m)
    R = range(N)
    for i in R:
        for j in R:
            if i != j:
                yield (i,j)
                ind += 1


if __name__ == "__main__":
    ## from pprint import pprint
    ## n = 7
    ## cmat = scipy.ones((n,n))
    ## for i, j in ((1,4),(1,7),
    ##              (2,6),(2,7),
    ##              (3,6),(3,7),
    ##              (4,5),
    ##              (5,6),):
    ##     cmat[i-1,j-1] = 0
    ##     cmat[j-1,i-1] = 0
    ## pprint(cmat)
    #pprint(list(nchoosem.dists_by_maxsize(n,4)))


## for x in nondiag_indices(scipy.zeros([6,6])):
##     print x
## sys.exit()

## m = RateModelGE(4)
## s = set()
## for i, d in m.enumerate_dists():
##     for as in iter_ancsplits(d):
##         if as.descdists not in s:
##             s.add(as.descdists)
##         else:
##             print "!", as.descdists

    import time
    t = time.time()
    m1 = RateModelGE(8)
    print time.time() - t, "to create v1"
    t = time.time()
    m2 = RateModelGE2(8)
    print time.time() - t, "to create v2"
    from pprint import pprint
    ## print list(m1.iter_dist_splits((0,1,1)))
    ## print list(m2.iter_dist_splits((1,2)))
    ## pprint(m1.Q)
    ## pprint(m2.Q)
    t = time.time()
    pprint(m1.P(0, 1.33))
    print time.time() - t
    t = time.time()
    pprint(m2.P(0, 1.33))
    print time.time() - t
    ## import scipy.sparse
    ## print m.Q[0].shape
    ## coo = scipy.sparse.coo_matrix(m.Q[0])
    ## print coo
    ## sys.exit()
    ## t = time.time()
    ## pprint(rates.Q2P(m.Q[0], 1.0))
    ## print time.time() - t, "for Q2P"
    ## t = time.time()
    ## pprint(rates.q2p_test(m.Q[0], 1.0))
    ## print time.time() - t, "for q2p_test"
