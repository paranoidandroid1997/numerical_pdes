# Applies boundary conditions to a function
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