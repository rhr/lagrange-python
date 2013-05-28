import sys, random, time
import scipy
from scipy import sum
import scipy.linalg
array = scipy.array
transpose = scipy.transpose
exp = scipy.exp
identity = scipy.identity
dot = scipy.dot
product = scipy.product
log = scipy.log
zeros = scipy.zeros
linalg = scipy.linalg

#from LinearAlgebra import eigenvectors, eigenvalues, inverse
#import phylo, newick

rand = random.Random()
uniform = rand.random

def Q2P_old(Q, t):
    """
    Q: a square rank-2 Numeric array representing the instantaneous rate
        matrix
    
    t: time interval (branch length)

    returns: P, the probability matrix of observing change in time t
        given values in Q

    ***********************************************************
    * This function computes the transformation described in: *
    * Pagel (1994). Proc. R. Soc. Lond. B 255: 37-45.         *
    * Please cite accordingly!                                *
    ***********************************************************

                    -1
    P(t) = Cexp(Dt)C  ,

    where C contains the eigenvectors of Q and Dt is a diagonal matrix
    of the form:

        |exp(lambda1*t)        0       |
    D = |                              |,
        |      0        exp(lambda2 *t)|

    where lambda{1,2} are the eigenvalues of Q.
    """
    evals, evecs = linalg.eig(Q)
    #C = transpose(evecs) # need to transpose these for dot product
    C = evecs
    exp_D = exp(evals*t)*identity(Q.shape[0])
    C_inv = linalg.inv(C)
    P = dot(dot(C, exp_D), C_inv)
    Psum = scipy.sum(P,1) # sum across rows, and
    P = P/Psum # divide by sum so that all elements are between 0 and 1
    return P

def Q2P(Q, t):
    p = linalg.expm(Q*t)
    return p/sum(p,1)

def binP(pi0, t):
    pi1 = 1.0 - pi0
    mu = 1.0/(2.0 * pi0 * pi1)
    eut = exp(-(mu * t))
    return array([
        [pi0 + (pi1 * eut), pi1 - (pi1 * eut)],
        [pi0 - (pi0 * eut), pi1 + (pi0 * eut)],
        ])


def Q2(pi0=0.5):
    """
    Rate matrix for character with two states, s0 and s1.  pi0 and pi1
    are the stationary frequencies of s0 and s1.
    """
    assert 0.0 < pi0 < 1.0
    pi1 = 1.0 - pi0
    Q = array([(-pi1, pi1),
               (pi0, -pi0)])
    Q *= 1.0/(2.0*pi0*pi1)
    return Q

def Q2Pdict(nodes, Q):
    d = {}
    for n in nodes:
        if n.parent:
            length = n.length
            d[length] = Q2P(Q, length)
    return d

def binPdict(nodes, pi0):
    d = {}
    for n in nodes:
        if n.parent:
            brlen = n.length
            d[brlen] = binP(pi0, brlen)
    return d

def fractionals(nodes, postorder, data, states, prior, Q2Pdict):
    nstates = len(states)
    zeroed_array = zeros((nstates,), "fd").copy
    
    nodeFractionalProbs = {}

    # iterate over nodes in post-order traversal
    for nodeNum in postorder:

        # n is the current node
        n = nodes[nodeNum]

        if n.istip:
            node_conds = zeroed_array()
            state = data[n.nodeNum]
            try:
                state = int(state)
                node_conds[state] = 1.0
            except ValueError:
                if state == '?' or state == '-':
                    node_conds = node_conds + 1/float(nstates)
            
        else:
            node_conds = zeroed_array()

            for s0 in states:
                s0_cond = []
                for child in n.children():
                    Ps0 = Q2Pdict[child.length][s0]
                    cc = nodeFractionalProbs[child.nodeNum]
                    s0_cond.append(sum([ Ps0[s1]*cc[s1] for s1 in states ]))
                node_conds[s0] = product(s0_cond)
            node_conds = array(node_conds)*prior

        fracs = array(node_conds)/sum(node_conds)

        nodeFractionalProbs[nodeNum] = fracs

    return nodeFractionalProbs


def sample_ancstates(nodes, preorder, states, fractionals, Qmap):
    """
    Sample ancestral states from their conditional probabilities.
    Return a mapping of nodeNum -> ancstate
    """
    ancstates = {}
    for n in [ nodes[i] for i in preorder ]: #if not nodes[i].istip ]:
        nodefracs = fractionals[n.nodeNum]

        if n.parent:
            P = Qmap[n.length]
            ancst = ancstates[n.parent.nodeNum]
            newstate_Prow = P[ancst]
            nodefracs = nodefracs * newstate_Prow
            nodefracs /= sum(nodefracs)

        rv = uniform()
        v = 0.0
        for state, frac in zip(states, nodefracs):
            v += frac
            if rv < v:
                break
        ancstates[n.nodeNum] = state

    return ancstates

if __name__ == "__main__":
    Q = Q2(0.5)
    ## print Q2P(Q, 0.5)
    ## print
    print Q2P(Q2(3.16060279e-01), 0.5)
    print
    print q2p_test(Q2(3.16060279e-01), 0.5)
    #print Q2P(Q, 1.0)
##     from pprint import pprint

##     t = newick.parse("(A:0.1, (B:0.1, C:0.1):0.1);")
##     phylo.polarize(t)
##     nodes = [ n for n in t.descendants() ]
##     for i, n in enumerate(nodes):
##         n.nodeNum = i
##     postorder = [ n.nodeNum for n in t.descendants(phylo.POSTORDER) ]
##     states = (0, 1)
##     prior = (1.0, 1.0)
##     data = {"A": 0, "B": 0, "C": 1}
##     for n in t.leaves():
##         data[n.nodeNum] = data[n.label]

##     Q = Q2()
##     Qmap = Q2Pdict(nodes, Q)

##     tfracs = fractionals(nodes, postorder, data, states, prior, Qmap)
##     pprint(tfracs)

##     ancs = sample_ancstates(nodes, range(len(nodes)), states, tfracs, Qmap)
##     print ancs

