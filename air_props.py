import numpy as np
from CoolProp.CoolProp import PropsSI


def calculate(T, P):
    T = np.asarray(T, dtype=np.float64).ravel()
    P = np.asarray(P, dtype=np.float64).ravel()

    cp = PropsSI("Cpmass", "T", T, "P", P, "Air")
    rho = PropsSI("Dmass", "T", T, "P", P, "Air")
    mu = PropsSI("V", "T", T, "P", P, "Air")
    k = PropsSI("L", "T", T, "P", P, "Air")

    # shape = (4, N)
    return np.stack((rho, cp, mu, k))