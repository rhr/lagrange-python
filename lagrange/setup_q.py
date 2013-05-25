def setup_Q(nperiods, dists, D, Dmask):
    t = 0.0
    ndists = len(dists)
    self.Q = scipy.zeros((nperiods, ndists, ndists))
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
