import matplotlib.pyplot as plt

from FluxEval import *
from applyBC import *

# Calculate u up to tmax given a method type
def evaluate(u, a, dt, dx, tmax, ibeg, iend, ngc, N, methodType, BCtype):
    t = 0

    # Initialize two variables to hold U^{n + 1} and U^{n}
    uNew = u.copy()

    while t < tmax:
        if methodType == 1:
            # Initialize flux
            Flux = np.zeros(N + (2 * ngc) + 1)

            # Calculate flux
            for i in range(ibeg, (iend + 2)):
                Flux[i] = FluxEval(
                    u[i - 2], u[i - 1], u[i], u[i + 1], a, dt, dx, methodType
                )

            for i in range(ibeg, (iend + 1)):
                uNew[i] = u[i] - ((dt / dx) * (Flux[i + 1] - Flux[i]))

        # Apply Boundary Conditions in place
        applyBC(uNew, ibeg, iend, ngc, BCtype)

        # Set new to old
        u = uNew.copy()

        # Increment time
        t += dt


def double_plot(x, y1, y2):
    fig, axs = plt.subplots(2,1)

    axs[0].plot(x, y1)
    axs[1].plot(x, y2)

    plt.show()


def runAndPlot(
    x, u1, u2, a, dt, dx, tmax1, tmax2, ibeg, iend, ngc, N, methodType, BCtype
):
    evaluate(u1, a, dt, dx, tmax1, ibeg, iend, ngc, N, methodType, BCtype)
    evaluate(u2, a, dt, dx, tmax2, ibeg, iend, ngc, N, methodType, BCtype)
    double_plot(x, u1, u2)
