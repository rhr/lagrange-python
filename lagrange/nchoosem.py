#!/usr/bin/env python
import sys, scipy

def comb(m, n):
    """m, n -> number of combinations of m items, n at a time.

    m >= n >= 0 required.
    """

    if not m >= n >= 0:
        raise ValueError("m >= n >= 0 required: " + `m, n`)
    if n > (m >> 1):
        n = m-n
    if n == 0:
        return 1
    result = long(m)
    i = 2
    m, n = m-1, n-1
    while n:
        # assert (result * m) % i == 0
        result = result * m / i
        i = i+1
        n = n-1
        m = m-1
    return result

def comb_at_index(m, n, i):
    """m, n, i -> i'th combination of m items taken n at a time.

    m >= n >= 1 and 0 <= i < comb(m, n) required.

    Return the i'th combination in lexicographic order, as a list
    of n elements taken from range(m).
    The index (i) is 0-based.

    Example:
    >>> for i in range(6):
    ...    print comb_at_index(4, 2, i)
    [0, 1]
    [0, 2]
    [0, 3]
    [1, 2]
    [1, 3]
    [2, 3]
    """

    if not m >= n >= 1:
        raise ValueError("m >= n >= 1 required: " + `m, n`)
    c = long(comb(m, n))
    if not 0 <= i < c:
        raise ValueError("0 <= i < comb(m,n) required: " + `i, c`)
    result = []
    # have c == comb(m, n), want comb(m-1,n-1)
    c = c * n / m
    # invariant: c == comb(m-1, n-1)
    for element in xrange(m):
        if i < c:
            # take this element, and n-1 from the remaining
            result.append(element)
            n = n-1
            if n == 0:
                break
            # have c == comb(m-1,n), want comb(m-2,n-1)
            c = c * n / (m-1)
        else:
            # skip this element, and take all from the remaining
            i = i-c
            # have c == comb(m-1,n-1), want comb(m-2,n-1)
            c = c * (m-n) / (m-1)
        m = m-1
    assert i == 0
    return result

def iterate(M, N):
    if not M >= N >= 1:
        raise ValueError("m >= n >= 1 required: " + `M, N`)
    ncombs = long(comb(M, N))
    for x in xrange(ncombs):
        i = x; n = N; m = M
        c = ncombs * n / m
        result = []
        for element in xrange(m):
            if i < c:
                # take this element, and n-1 from the remaining
                result.append(element)
                n = n-1
                if n == 0:
                    break
                # have c == comb(m-1,n), want comb(m-2,n-1)
                c = c * n / (m-1)
            else:
                # skip this element, and take all from the remaining
                i = i-c
                # have c == comb(m-1,n-1), want comb(m-2,n-1)
                c = c * (m-n) / (m-1)
            m = m-1
        assert i == 0
        yield tuple(result)

## for x in iterate(7,2):
##     print x
## sys.exit()

def dists_by_maxsize(nareas, maxsize):
    v = [0]*nareas
    for i in range(nareas):
        x = v[:]; x[i] = 1
        yield tuple(x)
    n = 2
    while n <= maxsize:
        for indices in iterate(nareas, n):
            x = v[:]
            for i in indices:
                x[i] = 1
            yield tuple(x)
        n += 1

def dists_by_maxsize_idx(nareas, maxsize):
    yield tuple()
    for i in range(nareas):
        yield (i,)
    n = 2
    while n <= maxsize:
        for indices in iterate(nareas, n):
            yield indices
        n += 1

## for x in dists_by_maxsize_idx(7, 2):
##     print x
## sys.exit()

def idx2bitvect(indices, M):
    v = [0]*M
    for i in indices: v[i] = 1
    return tuple(v)

def iterate_all(m):
    for n in range(1, m+1):
        it = iterate(m, n)
        for x in it:
            yield x

def iterate_all_bv(m):
    it = iterate_all(m)
    for x in it:
        yield idx2bitvect(x, m)

def iterate_all_bv2(m):
    yield tuple([0]*m)
    it = iterate_all(m)
    for x in it:
        yield idx2bitvect(x, m)

def iterate_all_idx(m):
    yield tuple()
    for x in iterate_all(m):
        yield x

def main():
    m, n = 4, 3
    
    #print "%d choose %d:" % (m, n),

    ## ncombs = comb(m, n)
    ## print "%s combinations" % ncombs

    ## for i, x in enumerate(iterate_all_bv2(m)):
    ##     print i+1, x

    for i, x in enumerate(iterate_all_idx(m)):
        print i, x
    
##     for i in range(ncombs):
##         cmb = comb_at_index(m, n, i)
##         print cmb

def subsets(s, size):
    m = len(s)
    for i in range(comb(m, size)):
        indices = comb_at_index(m, size, i)
        yield [ s[j] for j in indices ]

if __name__ == "__main__":
    ## s = "abcdef"
    ## print list(subsets(s, 2))
    print len(list(dists_by_maxsize_idx(12,4)))
