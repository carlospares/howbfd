def upwind(uL, uR, crit):
    return uL if crit>=0 else uR