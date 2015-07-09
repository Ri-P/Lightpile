# -*- coding: utf-8 -*-
"""
Calculate basic properties of plane waves.
"""
from __future__ import unicode_literals
from __future__ import print_function
import numpy as np

__author__ = "Richard Pfeifer"


def Sz(u, d, kz, k_vac, N, pol):
    """
    Calculate z-component of Poyntingvector.

    Unit: Sz * 2 * c0 * mu0 in SI-units with appropriate field amplitudes.
    """

    if pol == "TE":
        Sz_ar = (((kz.real *
                  (np.conjugate(u)*u - np.conjugate(d)*d).real) -
                  ((2. * kz).imag * (u*np.conjugate(d)).imag)) /
                 k_vac).astype(np.float32)  # 1/(2 c0 mu0) missing
    elif pol == "TM":
        nkz = kz * np.conjugate(N) / N
        Sz_ar = ((nkz.real *
                  (np.conjugate(u)*u - np.conjugate(d)*d).real -
                  2 * nkz.imag * (u*np.conjugate(d)).imag) /
                 k_vac).astype(np.float32)  # 1/(2 c0 mu0) missing
    else:
        raise ValueError("Invalid polarization %s" % pol)

    return Sz_ar


def kz(k_vac, N, q):
    """
    Calculate kz-component of k-vector.

    The time harmonic convention used forces kz.imag <= 0
    """
    kz = np.sqrt((k_vac * N)**2 - q**2)
    return np.where(kz.imag <= 0, kz, -kz)


def w(N, u):
    """
    Calculate normed transversal component of k-vector.
    """
    w = np.sqrt(N**2 - u**2)
    return np.where(w.imag <= 0, w, -w)


def field_squared(u, d):
    """
    Calculate the squared amplitude of the total electric field.
    """

    return abs((u + d).conjugate() * (u + d))
