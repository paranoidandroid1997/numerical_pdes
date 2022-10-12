import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from applyBC import *

# Set grid resolution
N = 64  # input("Grid resolution N = ")

# Set x bounds
xa = 0
xb = 1

# Set change x
dx = (xb - xa) / N

# Set boundaries
ngc = 2
ibeg = ngc
iend = (N + ngc) - 1

# Set values for the left and right values of the square wave
uL = -1
uR = 2

# Intialize x array
x = np.zeros(N + (2 * ngc))
x[ibeg : (iend + 1)] = np.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

# Initialize boundaries of x
for i in range(0, ngc):
    x[ibeg - 1 - i] = x[ibeg] - (i + 1) * dx
    x[i + N + ngc] = x[N + ngc - 1] + (i + 1) * dx

# Initialize u and uInit
u = np.zeros(N + (2 * ngc))
uInit = np.zeros(N + (2 * ngc))

# Initialize flux
flux = np.zeros(N + (2 * ngc) + 1)

# IC setup
shockLoc = 0.5
travelDist = 0.1

# Setup square wave
for i in range(ibeg, iend + 1):
    if x[i] > xa and x[i] < shockLoc * (xb - xa):
        u[i] = uL
    elif x[i] > shockLoc * (xb - xa):
        u[i] = uR

for i in range(ibeg, iend):
    if x[i] > xa and x[i] <= (shockLoc + travelDist) * (xb - xa):
        uInit[i] = uL
    elif x[i] > (shockLoc + travelDist) * (xb - xa):
        uInit[i] = uR

# Get boundary condition type
BCtype = 2

# Apply BC in place
applyBC(u, ibeg, iend, ngc, BCtype)
applyBC(uInit, ibeg, iend, ngc, BCtype)

# Set cfl condition
cfl = 0.9  # input("Cfl = ")

# Set timestep

dt = cfl * dx / abs(np.max(u))
Ncycle = 1
t = 0

tmax = 1  # input("tmax=")

# Plot the initial state of u
# plt.plot(x[ibeg:(iend + 1)], u[ibeg:(iend+1)])
# plt.show()


for i in range(ibeg, iend + 2):
    # Compute shock speed for Burgers' eqn
    s = 0.5 * (u[i] + u[i - 1])

    # Upwind Flux
    if (s < 0) and (u[i] < 0):
        flux[i] = 0.5 * u[i] ** 2
    elif (s > 0) and (u[i - 1] > 0):
        flux[i] = 0.5 * u[i - 1] ** 2
    elif (u[i - 1] < 0) and (u[i] > 0):
        flux[i] = 0


def animateConservative(i):
    global u

    uNew = np.zeros(N + (2 * ngc))

    for i in range(ibeg, iend + 2):
        # Compute shock speed for Burgers' eqn
        s = 0.5 * (u[i] + u[i - 1])

        # Upwind Flux
        if (s < 0) and (u[i] < 0):
            flux[i] = 0.5 * u[i] ** 2
        elif (s > 0) and (u[i - 1] > 0):
            flux[i] = 0.5 * u[i - 1] ** 2
        elif (u[i - 1] < 0) and (u[i] > 0):
            flux[i] = 0

    for i in range(ibeg, iend + 1):
        uNew[i] = u[i] - dt / dx * (flux[i + 1] - flux[i])

    applyBC(uNew, ibeg, iend, ngc, BCtype)

    u = uNew

    line.set_data(x, u)

    return (line,)


def animateNonConservative(i):
    global u

    uNew = np.zeros(N + (2 * ngc))

    for i in range(ibeg, iend + 1):
        uNew[i] = u[i] - u[i] * dt / dx * (u[i] - u[i - 1])

    applyBC(uNew, ibeg, iend, ngc, BCtype)

    u = uNew

    line.set_data(x, u)

    return (line,)


fig, ax = plt.subplots()
ax.set_xlim(-0.2, 1.2)
ax.set_ylim(-2, 2)
ax.grid()
(line,) = ax.plot([], [], lw=2, marker="o", markerfacecolor="none", markeredgecolor="r")

u0 = u.copy()
anim = FuncAnimation(fig, animateConservative, frames=200, interval=20, blit=True)
anim.save("./animations/Conservative.gif", fps=60, dpi=200)

u = u0.copy()
anim = FuncAnimation(fig, animateNonConservative, frames=120, interval=20, blit=True)
anim.save("./animations/NonConservative.gif", fps=60, dpi=200)
