import numpy as np
import matplotlib.pyplot as plt

from applyIC import *
from applyBC import *
from FluxEval import *
from runAndPlot import *


def quick_plot(x, y, xlims, ylims, title="placeholder", show=True, save=False):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.grid(alpha=0.33)
    ax.set_title(title)

    ax.plot(
        x,
        y,
        color="slateblue",
        marker="o",
        ms=5,
        markerfacecolor="None",
        markeredgecolor="lime",
        markeredgewidth=1,
    )

    if save:
        fig.save(f"{title}.png")
    if show:
        plt.show()


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

    # Initial conditions plot
    # quick_plot(x, u, xlims, ylims)

    print("[1]Upwind, [2]LW, [3]Fromm, [4]BW")
    print("[5]minmod, [6]superbee, [7]MC, [8]VanLeer")
    print("[9]LF")

    methodType = 1  # int(input("Method type [1-9] = "))
    runAndPlot(
        x, u1.copy(), u2.copy(), uInit1, uInit2, a, dt, dx, tmax1, tmax2, ibeg, iend, ngc, N, methodType, BCtype
    )

    methodType = 2  # int(input("Method type [1-9] = "))
    runAndPlot(
        x, u1.copy(), u2.copy(), uInit1, uInit2, a, dt, dx, tmax1, tmax2, ibeg, iend, ngc, N, methodType, BCtype
    )

    # while t < tmax:

    #     # Calculate flux
    #     #for i in range(ibeg, (iend + 2)):
    #         #Flux[i] = FluxEval(u[i-2], u[i-1], u[i], u[i+1], a, dt, dx, methodType)

    #     # Calculate uNew
    #     for i in range(ibeg, (iend + 1)):
    #         if (methodType == 1):
    #             # Upwind
    #             uNew[i] = u[i] - ((dt/dx) * (Flux[i + 1] - Flux[i]))
    #         elif (methodType == 2):
    #             # LW
    #             uNew[i] = u[i] - (a/2)*(dt/dx)*(u[i + 1] - u[i - 1]) + (1/2)*(a*dt/dx)**2*(u[i + 1] - 2*u[i] + u[i - 1])
    #         elif (methodType == 9):
    #             # LF
    #             uNew[i] = 0.5 * (u[i + 1] + u[i - 1]) - ((a/2) * (dt/dx) * (u[i + 1] - u[i - 1]))
    #         elif (methodType == 4):
    #             # BW
    #             uNew[i] = u[i] - (a/2)*(dt/dx)*(3*u[i] - 4*u[i-1] + u[i-2]) + (1/2)*(a*dt/dx)**2*(u[i] - 2*u[i-1] + u[i-2])
    #         elif (methodType == 3):
    #             # Fromm
    #             uNew[i] = u[i] - a * (dt/dx)*(u[i] - u[i - 1]) - (a/4)*(dt/dx)*(1 - a*(dt/dx))*(u[i+1] - u[i]) + (a/4)*(dt/dx)*(1 - a*(dt/dx))*(u[i-1] - u[i-2])

    #     # Apply Boundary Conditions in place
    #     applyBC(uNew, ibeg, iend, ngc, BCtype)

    #     # Set new to old
    #     u = uNew.copy()

    #     # Increment time
    #     t += dt

    #     quick_plot(x, u, xlims, ylims)

    # fig, ax = plt.subplots(figsize = (12, 6))

    # ax.set_xlim(xlims)
    # ax.set_ylim(ylims)
    # ax.grid(alpha=0.33)
    # ax.set_title("placeholder")

    # ax.plot(x, u, color ='slateblue',marker='o', ms = 5, markerfacecolor = "None",
    #         markeredgecolor="lime", markeredgewidth=1,)

    # ax.plot(x,uInit, color = "red")

    # plt.show()
