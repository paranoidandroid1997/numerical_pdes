# Not entirely sure what this does yet
def minmod(a, b):
    if ((a*b) > 0):
        if (np.abs(a) < np.abs(b)):
            c = a
        else:
            c = b
    else:
        c = 0
    return c