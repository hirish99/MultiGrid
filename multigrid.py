from interpolation import generate_interpolation
from poisson import get_poisson_fd_3d



def jac_smoothing(phi, f, h, n,nsmooth):
    omega = 2/3
    delta_x = h
    delta_y = h 
    delta_z = h
    di = (2/delta_x**2) + (2/delta_y**2) + (2/delta_z**2) 
    A = get_poisson_fd_3d(n)
    uk = phi
    for i in range(nsmooth):    
        uk = uk + omega * 1/di * (f - A @ phi)
        r = f - A @ uk
    return uk


def V_cycle(phi, f, h):
    # Pre-smoothing
    phi = jac_smoothing(phi, f, h)