import numpy as np
from poisson import get_poisson_fd_3d
import scipy


def jac_smoothing(f, h, nsmooth, A, n):
    omega = 2/3
    di = h**2/6
    uk = np.zeros(n**3)
    r = f - A @ uk
    for i in range(nsmooth):    
        uk = uk + omega * di * r
        r = f - A @ uk
    return uk

def coarse_grid_correction(n, ):


def test_smoother():
    n = 4 #should be 1 less than a power of 2
    h = 1/(n+1)
    print(get_poisson_fd_3d(n).todense())
    A = get_poisson_fd_3d(n) * 1/h**2
    #k = np.linspace(0,1,n)
    #V_1D = np.sin((h*np.pi)*np.outer(k,k))
    #V_3D = np.sqrt(2*h) * np.sqrt(2*h) * np.sqrt(2*h) * scipy.sparse.kron(scipy.sparse.kron(V_1D,V_1D),V_1D)
    ue = np.random.rand(n**3)
    #ue = V @ ue
    b = A @ ue

    us = jac_smoothing(b, h, 20, A, n)
    err = np.linalg.norm(b - A@us)
    print(err)

    

test_smoother()


