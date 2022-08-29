import numpy as np

if __name__ == "__main__":

    # Get number of steps in space from user
    N = int(input("Grid resolution N = "))

    # Set space bounds
    xa = 0
    xb = 1

    # Calculate needed change of x for the input grid resolution and space bounds
    dx = (xb - xa)/N

    # TBD what exactly this value is
    ngc = 2

    # Set the bounds of the interval (not counting BCs)
    ibeg = (ngc + 1) - 1
    iend = (N + ngc) - 1

    # Initialize main interval of x
    x = np.zeros(N + (2*ngc))
    x[ibeg:(iend + 1)] = np.linspace(xa+0.5*dx, xb-0.5*dx, N)

    # Initialize boundaries of x
    for i in range(0, ngc):
        x[i] = x[ibeg] - (i + 1)*dx
        x[i + N + ngc] = x[N + ngc - 1] + (i + 1)*dx

    # Initialize u
    u = np.zeros(N + (2*ngc))

    # Initialize flux
    Flux = np.zeros(N + (2*ngc) + 1)

