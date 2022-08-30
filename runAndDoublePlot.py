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
            # Upwind

            # Initialize flux
            Flux = np.zeros(N + (2 * ngc) + 1)

            # Calculate flux
            for i in range(ibeg, (iend + 2)):
                Flux[i] = FluxEval(
                    u[i - 2], u[i - 1], u[i], u[i + 1], a, dt, dx, methodType
                )

            for i in range(ibeg, (iend + 1)):
                uNew[i] = u[i] - ((dt / dx) * (Flux[i + 1] - Flux[i]))
        if methodType == 2:
            # LF
            for i in range(ibeg, (iend + 1)):
                uNew[i] = 0.5 * (u[i + 1] + u[i - 1]) - (
                    (a / 2) * (dt / dx) * (u[i + 1] - u[i - 1])
                )
        if methodType == 3:
            # LW
            for i in range(ibeg, (iend + 1)):
                uNew[i] = (
                    u[i]
                    - (a / 2) * (dt / dx) * (u[i + 1] - u[i - 1])
                    + (1 / 2) * (a * dt / dx) ** 2 * (u[i + 1] - 2 * u[i] + u[i - 1])
                )
        if methodType == 4:
            # BW
            for i in range(ibeg, (iend + 1)):
                uNew[i] = (
                    u[i]
                    - (a / 2) * (dt / dx) * (3 * u[i] - 4 * u[i - 1] + u[i - 2])
                    + (1 / 2) * (a * dt / dx) ** 2 * (u[i] - 2 * u[i - 1] + u[i - 2])
                )
        if methodType == 5:
            # Fromms
            for i in range(ibeg, (iend + 1)):
                uNew[i] = (
                    u[i]
                    - a * (dt / dx) * (u[i] - u[i - 1])
                    - (a / 4) * (dt / dx) * (1 - a * (dt / dx)) * (u[i + 1] - u[i])
                    + (a / 4) * (dt / dx) * (1 - a * (dt / dx)) * (u[i - 1] - u[i - 2])
                )

        # Apply Boundary Conditions in place
        applyBC(uNew, ibeg, iend, ngc, BCtype)

        # Set new to old
        u[:] = uNew.copy()[:]

        # Increment time
        t += dt


def doublePlot(x, y1, y2, reference1, reference2, title="placehodler", save=False):
    fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    axs[0].set_title(title)

    axs[0].plot(
        x,
        y1,
        linewidth=1,
        color="slateblue",
        marker="o",
        ms=3,
        markerfacecolor="None",
        markeredgecolor="red",
        markeredgewidth=1,
    )

    axs[0].plot(x, reference1, color="black")

    axs[1].plot(
        x,
        y2,
        linewidth=1,
        color="slateblue",
        marker="o",
        ms=3,
        markerfacecolor="None",
        markeredgecolor="red",
        markeredgewidth=1,
    )

    axs[1].plot(x, reference2, color="black")

    if save:
        fig.savefig(f"./plots/{title}.png", dpi=200)
    plt.show()


def runAndDoublePlot(
    x,
    u1,
    u2,
    reference1,
    reference2,
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
):
    methods = ["Upwind", "LF", "LW", "BW", "Fromm"]
    evaluate(u1, a, dt, dx, tmax1, ibeg, iend, ngc, N, methodType, BCtype)
    evaluate(u2, a, dt, dx, tmax2, ibeg, iend, ngc, N, methodType, BCtype)
    doublePlot(x, u1, u2, reference1, reference2, methods[methodType - 1], save=True)
