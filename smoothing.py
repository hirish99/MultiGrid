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

def test_smoother():
    n = 8
    h = 1/(n+1)
    A = get_poisson_fd_3d(n) * 1/h**2
    #k = np.linspace(0,1,n)
    #V_1D = np.sin((h*np.pi)*np.outer(k,k))
    #V = np.sqrt(2*h) * np.sqrt(2*h) * np.sqrt(2*h) * scipy.sparse.kron(scipy.sparse.kron(V_1D,V_1D),V_1D)
    ue = np.random.rand(n**3)
    #ue = V @ ue
    b = A @ ue

    us = jac_smoothing(b, h, 3, A, n)
    err = np.linalg.norm(b - A@us)
    print(err)

    

test_smoother()


