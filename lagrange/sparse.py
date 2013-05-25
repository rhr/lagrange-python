import nchoosem
import scipy.sparse

def create_sparse_matrix(nareas):
    N = 2**nareas
    m = scipy.sparse.lil_matrix((N,N))
    for i in range(N):
        pass
