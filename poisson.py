import numpy as np
import scipy

def get_poisson_fd_3d(n):
    diagonals = [np.ones(n-1), -2*np.ones(n), np.ones(n-1)]
    offsets = np.array([-1, 0, 1])
    Ax = scipy.sparse.diags(diagonals, offsets, format='csc')

    A2D = scipy.sparse.kron(np.eye(n), Ax) + scipy.sparse.kron(Ax, np.eye(n))
    A3D = scipy.sparse.kron(np.eye(n), A2D)  + scipy.sparse.kron(A2D,np.eye(n))

    return A3D









