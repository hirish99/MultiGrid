import numpy as np

def generate_interpolation(nc):
    no = 2*nc

    P = np.zeros((no-1, nc-1))

    k = 1
    for i in range(nc-1):
        P[k-1, i] = 0.5
        P[k, i] = 1
        P[k+1, i] = 0.5
        k += 2
    return P