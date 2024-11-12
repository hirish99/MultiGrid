from interpolation import generate_interpolation
from smoothing import jac_smoothing






def V_cycle(phi, f, h):
    # Pre-smoothing
    phi = jac_smoothing(phi, f, h)