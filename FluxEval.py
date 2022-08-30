# Evaluates the flux for a given u and method type
def FluxEval(uLL, uL, uR, uRR, a, dt, dx, methodType):

    # Set some constants
    epsilon = 1.e-16
    ap = max(a, 0)
    am = min(a, 0)
    aAbs = ap - am
    amod = np.abs(a)
    Ca = amod * (dt/dx)

    if (a > 0):
        delU1 = uL - uLL
    else:
        delU1 = uRR - uR

    if (methodType <= 3):
        delU2 = uR - uL
        theta = 0
    elif (methodType > 3):
        delU2 = max(uR - uL, epsilon)
        theta = delU1/delU2

    # The different methods
    if (methodType == 1):
        # Upwind
        phi = 0
    elif (methodType == 2):
        # LW
        phi = 1
    elif (methodType == 3):
        # Fromm
        phi = 0.5 * (1.0 + theta)
    elif (methodType == 4):
        # BW
        phi = theta
    elif (methodType == 5):
        # minmod
        phi = minmod(1, theta)
    elif (methodType == 6):
        # superbee
        phi = max(0, max(min(1,2 * theta), min(2, theta)))
    elif (methodType == 7):
        # MC
        phi = max(0, min(min(0.5 * (1 + theta),2), 2 * theta))
    elif (methodType == 8):
        # vanLeer
        phi = (theta + np.abs(theta)) / (1 + np.abs(theta))
    elif (methodType == 9):
        # LF placeholder
        phi = 0

    if (methodType == 3):
        delta = 0.5 * (delU2 + delU1)
    else:
        delta = phi * delU2
    
    Flux = (am * uR) + (ap *uL) + (0.5 * amod * (1 - Ca) * delta)
    return Flux