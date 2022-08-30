import numpy as np


def applyIC(u, x, xa, xb, ibeg, iend, ICtype):
    # Initialize initial condititions based on user input
    if ICtype == 1:
        # Initialize square wave initial conditions
        for i in range(ibeg, (iend + 1)):
            if (x[i] > xa) and (x[i] <= (0.25 * (xb - xa))):
                u[i] = 0
            elif (x[i] > (0.25 * (xb - xa))) and (x[i] <= (0.5 * (xb - xa))):
                u[i] = -1
            elif (x[i] > (0.5 * (xb - xa))) and (x[i] <= (0.75 * (xb - xa))):
                u[i] = 1
            else:
                u[i] = 0
    elif ICtype == 2:
        # Initialize sin wave initial conditions
        u[ibeg : (iend + 1)] = np.sin(2.0 * np.pi * x[ibeg : (iend + 1)])
    elif ICtype == 3:
        # Initialize single shock square wave initial condition
        for i in range(ibeg, (iend + 1)):
            travelDist = 0.3
            if (x[i] > xa) and (x[i] <= ((0.5) * (xb - xa))):
                u[i] = 1
            elif x[i] > ((0.5) * (xb - xa)):
                u[i] = -1
    else:
        print("That is not a valid input")
