import numpy as np
import matplotlib.pyplot as plt

from applyIC import *
from applyBC import *
from FluxEval import *
from runAndDoublePlot import *
from runAndSinglePlot import *

if __name__ == "__main__":

    # Get number of steps in space from user
    N = 64  # int(input("Grid resolution N = "))

    # Set space bounds
    xa = 0
    xb = 1

    # Calculate needed change of x for the input grid resolution and space bounds
    dx = (xb - xa) / N

    # TBD what exactly this value is
    ngc = 2

    # Set the bounds of the interval (not counting BCs)
    ibeg = (ngc + 1) - 1
    iend = (N + ngc) - 1

    # Initialize main interval of x
    x = np.zeros(N + (2 * ngc))
    x[ibeg : (iend + 1)] = np.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

    # Initialize boundaries of x
    for i in range(0, ngc):
        x[ibeg - 1 - i] = x[ibeg] - (i + 1) * dx
        x[i + N + ngc] = x[N + ngc - 1] + (i + 1) * dx

    # Initialize u1 and u2
    u1 = np.zeros(N + (2 * ngc))
    u2 = np.zeros(N + (2 * ngc))
    uInit1 = np.zeros(N + (2 * ngc))
    uInit2 = np.zeros(N + (2 * ngc))

    # Get type of initial condition from the user
    ICtype1 = (
        2  # int(input("What type of IC for u1? [1:square; 2:sine, 3:single shock] = "))
    )
    ICtype2 = (
        3  # int(input("What type of IC for u2? [1:square; 2:sine, 3:single shock] = "))
    )

    # Apply initial conditions in place
    applyIC(u1, x, xa, xb, ibeg, iend, ICtype1)
    applyIC(u2, x, xa, xb, ibeg, iend, ICtype2)

    # Some values for plotting
    xlims = ((xa - 0.25), (xb + 0.25))
    ylims = (np.min(u1) * 1.5, np.max(u2) * 1.5)

    # Get type of boundary condition from the user
    BCtype = 1  # int(input("BC type [1:periodic, 2:outflow] = "))

    # Apply boundary conditions in place
    applyBC(u1, ibeg, iend, ngc, BCtype)
    applyBC(u2, ibeg, iend, ngc, BCtype)

    # Make a reference u0
    u10 = u1.copy()
    u20 = u2.copy()

    # Setup a uInit for later plotting
    if ICtype1 == 3:
        travelDist = 0.3
        for i in range(ibeg, (iend + 1)):
            if (x[i] > xa) and (x[i] <= ((0.5 + travelDist) * (xb - xa))):
                uInit1[i] = 1
            elif x[i] > ((0.5 + travelDist) * (xb - xa)):
                uInit1[i] = -1
        applyBC(uInit1, ibeg, iend, ngc, BCtype)
    else:
        uInit1 = u1.copy()

    # Setup a uInit for later plotting
    if ICtype2 == 3:
        travelDist = 0.3
        for i in range(ibeg, (iend + 1)):
            if (x[i] > xa) and (x[i] <= ((0.5 + travelDist) * (xb - xa))):
                uInit2[i] = 1
            elif x[i] > ((0.5 + travelDist) * (xb - xa)):
                uInit2[i] = -1
        applyBC(uInit2, ibeg, iend, ngc, BCtype)
    else:
        uInit2 = u2.copy()

    # Get advection velocity and cfl from the user

    a = 1  # float(input("Advection velocity a = "))
    cfl = 0.9  # float(input("CFL = "))

    # Calculate dt
    dt = cfl * (dx / np.abs(a))

    # Choose number of cycles
    Ncycle = 1

    # Set tmax value
    if ICtype1 == 3:
        tmax1 = travelDist / a
    else:
        tmax1 = Ncycle * ((xb - xa) / np.abs(a))

    # Set tmax value
    if ICtype2 == 3:
        tmax2 = travelDist / a
    else:
        tmax2 = Ncycle * ((xb - xa) / np.abs(a))

    for i in range(1, 6):
        methodType = i  # int(input("Method type [1-9] = "))
        runAndDoublePlot(
            x,
            u1.copy(),
            u2.copy(),
            uInit1,
            uInit2,
            a,
            dt,
            dx,
            tmax1,
            tmax2,
            ibeg,
            iend,
            ngc,
            N,
            methodType,
            BCtype,
        )
