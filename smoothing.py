import numpy as np
from poisson import get_poisson_fd_3d

def jac_smoothing(f, h, nsmooth, A):
    omega = 2/3
    delta_x = h
    delta_y = h
    delta_z = h
    di = (2/delta_x**2) + (2/delta_y**2) + (2/delta_z**2)
    uk = 0
    r = f - A @ uk
    for i in range(nsmooth):    
        uk = uk + omega * 1/di * r
        r = f - A @ uk
    return uk

def test_smoother():
    n = 4
    A = get_poisson_fd_3d(n)
    print(A.shape)
    #ue = np.random.rand(n)
    #b = A @ ue

test_smoother()


