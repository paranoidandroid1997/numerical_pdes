import numpy as np
import matplotlib.pyplot as plt

def quick_plot(x, y, xlims, ylims, title = "placeholder", show=True, save=False):
    fig, ax = plt.subplots(figsize = (12, 6))

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.grid(alpha=0.33)
    ax.set_title(title)

    ax.plot(x, y, color ='slateblue',marker='o', ms = 5, markerfacecolor = "None",
            markeredgecolor="lime", markeredgewidth=1,)


    if (save):
        fig.save(f"{title}.png")
    if (show):
        plt.show()

def applyBC(u, ibeg, iend, ngc, BCtype):
    if (BCtype == 1):
        # Periodic boundary conditions
        u[ibeg - ngc] = u[iend - 1]
        u[ibeg - ngc + 1] = u[iend]
        u[iend + 1] = u[ibeg]
        u[iend + 2] = u[ibeg + 1]
    elif (BCtype == 2):
        # Outflow bounddary conditions
        u[1:ngc] = u[ibeg]
        u[(iend + 1):(iend + ngc)] = u[iend]

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
    ICtype = int(input("What type of IC? [1:square; 2:sine, 3:single shock = "))

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
            if ((x[i] > xa and x[i] <= (0.5 * (xb - xa)))):
                u[i] = 1
            elif (x[i] > 0.5 * (xb - xa)):
                u[i] = -1
    else:
        print("That is not a valid input")

    # Some values for plotting
    xlims = ((xa - 0.25), (xb + 0.25))
    ylims = (np.min(u) * 1.2, np.max(u) * 1.2)

    # Debug plot
    quick_plot(x, u, xlims, ylims)

    # Finish setting up initial conditions
    if (ICtype == 3):
        travelDist = 0.3
        for i in range(ibeg, (iend + 1)):
            if ((x[i] > xa) and (x[i] <= ((0.5 + travelDist) * (xb - xa)))):
                    uInit[i] = 1
            elif (x[i] > ((0.5 + travelDist) * (xb - xa))):
                uInit[i] = -1
    else:
        uInit = u


    # Debug plot
    # quick_plot(x, uInit, xlims, ylims)


    # Get type of boundary condition from the user
    BCtype = int(input("BC type [1:periodic, 2:outflow] = ")) 

    # Apply boundary conditions in place
    applyBC(u, ibeg, iend, ngc, BCtype)

    print(x)

    # Final plot
    quick_plot(x, u, xlims, ylims)
