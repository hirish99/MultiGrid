import numpy as np
from poisson import get_poisson_fd_3d

def jac_smoothing(f, h, nsmooth, A):
    omega = 2/3
    di = 6/h**2
    uk = 0
    r = f - A @ uk
    for i in range(nsmooth):    
        uk = uk + omega * 1/di * r
        r = f - A @ uk
    return uk

def test_smoother():
    n = 4
    h = 1/(n+1)
    A = get_poisson_fd_3d(n) * 1/h**2
    #print(A.shape)
    ue = np.random.rand(n)
    b = A @ ue

    

test_smoother()


