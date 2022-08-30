import numpy as np
import matplotlib.pyplot as plt

from applyBC import *
from minmod import *
from FluxEval import *

def quick_plot(x, y, xlims, ylims, title = "placeholder", show=True, save=False):
    fig, ax = plt.subplots(figsize = (12, 6))

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.grid(alpha=0.33)
    ax.set_title(title)

    ax.plot(x, y, color ='slateblue',marker='o', ms = 5, markerfacecolor = "None",
            markeredgecolor="lime", markeredgewidth=1,)


    ax.plot(x,uInit, color = "red")

    if (save):
        fig.save(f"{title}.png")
    if (show):
        plt.show()

def double_plot(x, y1, y2,titles=("placeholder","placeholder"),
                show=True, save=False):
    fig, axs = plt.subplots(2, 1)

    axs[0].set_xlim(xlims)
    axs[0].set_ylim(ylims)
    axs[0].grid(alpha=0.33)
    axs[0].set_title(titles[0])

    axs[0].set_xlim(xlims)
    axs[0].set_ylim(ylims)
    axs[0].grid(alpha=0.33)
    axs[0].set_title(titles[1])

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
        x[ibeg - 1 - i] = x[ibeg] - (i + 1)*dx
        x[i + N + ngc] = x[N + ngc - 1] + (i + 1)*dx

    # Initialize u
    u = np.zeros(N + (2*ngc))
    uInit = np.zeros(N + (2*ngc))

    # Initialize flux
    Flux = np.zeros(N + (2*ngc) + 1)

    # Get type of initial condition from the user
    ICtype = int(input("What type of IC? [1:square; 2:sine, 3:single shock] = "))

    # Initialize initial condititions based on user input
    if (ICtype == 1):
        # Initialize square wave initial conditions
        for i in range(ibeg, (iend + 1)):
            if ((x[i] > xa) and (x[i] <= (0.25 * (xb - xa)))):
                u[i] = 0
            elif ((x[i] > (0.25 * (xb - xa))) and (x[i] <= (0.5 * (xb - xa)))):
                u[i] = -1
            elif ((x[i] > (0.5 * (xb - xa))) and (x[i] <= (0.75 * (xb - xa)))):
                u[i] = 1
            else:
                u[i] = 0
    elif (ICtype == 2):
        # Initialize sin wave initial conditions
        u[ibeg:(iend + 1)] = np.sin(2.0 * np.pi *x[ibeg:(iend + 1)])
    elif (ICtype == 3):
        #Initialize single shock square wave initial condition
        for i in range(ibeg, (iend + 1)):
            travelDist = 0.3
            if ((x[i] > xa) and (x[i] <= ((0.5) * (xb - xa)))):
                u[i] = 1
            elif (x[i] > ((0.5) * (xb - xa))):
                u[i] = -1
    else:
        print("That is not a valid input")

    # Some values for plotting
    xlims = ((xa - 0.25), (xb + 0.25))
    ylims = (np.min(u) * 1.5, np.max(u) * 1.5)

    # Debug plot
    #quick_plot(x, u, xlims, ylims)

    # Finish setting up initial conditions
    if (ICtype == 3):
        travelDist = 0.3
        for i in range(ibeg, (iend + 1)):
            if ((x[i] > xa) and (x[i] <= ((0.5 + travelDist) * (xb - xa)))):
                    uInit[i] = 1
            elif (x[i] > ((0.5 + travelDist) * (xb - xa))):
                uInit[i] = -1
    else:
        uInit = u.copy()

    # Debug plot
    # quick_plot(x, uInit, xlims, ylims)


    # Get type of boundary condition from the user
    BCtype = int(input("BC type [1:periodic, 2:outflow] = ")) 

    # Apply boundary conditions in place
    applyBC(u, ibeg, iend, ngc, BCtype)

    # Initialize two variables to hold U^{n + 1} and U^{n}
    uNew = u.copy()
    uOld = u.copy()

    # Get advection velocity and cfl from the user
    a = float(input("Advection velocity a = "))
    cfl = float(input("CFL = "))

    # Calculate dt
    dt = cfl * (dx/np.abs(a))

    # Choose number of cycles
    Ncycle = 1

    # Set initial t value and find a tmax
    t = 0
    if (ICtype == 3):
        tmax = travelDist/a
    else:
        tmax = Ncycle * ((xb - xa)/np.abs(a))

    # Initial conditions plot
    #quick_plot(x, u, xlims, ylims)









    # print("[1]Upwind, [2]LW, [3]Fromm, [4]BW")
    # print("[5]minmod, [6]superbee, [7]MC, [8]VanLeer")
    # print("[9]LF")
    # methodType = int(input("Method type [1-9] = "))

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

