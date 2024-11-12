import numpy as np
import scipy

def get_poisson_fd_3d(n):
    diagonals = [-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)]
    offsets = np.array([-1, 0, 1])
    Ax = scipy.sparse.diags(diagonals, offsets, format='csc')

    A2D = scipy.sparse.kron(np.eye(n), Ax) + scipy.sparse.kron(Ax, np.eye(n))
    A3D = scipy.sparse.kron(np.eye(n), A2D)  + scipy.sparse.kron(Ax, np.eye(n**2)) 

    return A3D

def get_poisson_fd_3d_jacobi(n):
    diagonals = [np.ones(n-1), 0*np.ones(n), np.ones(n-1)]
    offsets = np.array([-1, 0, 1])
    Ax = scipy.sparse.diags(diagonals, offsets, format='csc')

    A2D = scipy.sparse.kron(np.eye(n), Ax) + scipy.sparse.kron(Ax, np.eye(n))
    A3D = scipy.sparse.kron(np.eye(n), A2D)  + scipy.sparse.kron(Ax, np.eye(2*n)) 

    return A3D

#A3D = get_poisson_fd_3d_jacobi(2)
#print(A3D.todense())








