# -*- coding: utf-8 -*-
"""
Calculate basic properties of a dipole emitter.
"""
from __future__ import unicode_literals
from __future__ import print_function
import numpy as np

__author__ = "Richard Pfeifer"


norm = np.sqrt(3. / (8. * np.pi))


def dipole_amplitude_vw_d2v(dp_ori, dp_pol, v, w, n_dp):
    """
    Calculate a(v,w) in units of d^2v

    v =  q / k_vac
    w = kz / k_vac  = sqrt(n_dp^2 - v^2)   for v^2 <= n_dp^2 and n_dp.imag=0

    """
    if dp_ori == "v" and dp_pol == "TE":
        return (0, 0)
    elif dp_ori == "v" and dp_pol == "TM":
        a = - v * 1j * norm / (n_dp * w)
        return (a, a)

    else:
        raise ValueError("Invalid polarization '%s' or orientation '%s'" %
                         (dp_pol, dp_ori))


def dipole_amp_qkz(dp_ori, dp_pol, k_dp, k_z, q, phi=0.):
    """
    Get field amplitudes of a free dipole in medium dp with normed wavevector
        k_dp = k_vac * n_dp
    and direction k_dp = (q, k_z) where q is the projection of k_dp onto the
    kx-ky plane.
    k_z is assumed positive and the angle phi = angle(x-axis, q) = 0.

    Returned are the amplitudes u_0, d_0 in a plane above, below the emitter
    as a pair (u_0, d_0).

    (u_0, d_0) are normed field amplitudes.
    The actual field amplitudes in SI units u_SI, d_SI are

               sqrt( 8 pi)                      k_vac^2 p
       u_SI = -----------  A u_0    with A = ------------------
               sqrt(  3  )                    4 pi epsilon_vac

    Given the  amplitudes u_0 and d_0 integrated over all directions the normed
    emitted power of a vertical dipole in a homogeneous medium with purely real
    refractive index n is calculated to be

        P_hom_normed = Int(Sz(u_0, d_0) d^2q ) = n

    The unnormed power of this dipole in SI-units is
                          1              8 pi
        P_hom_SI = --------------- |A|^2 -----  P_hom_normed
                    2 c_vac mu_vac         3

    """
    if dp_ori == "v" and dp_pol == "TE":
        return (0, 0)
    elif dp_ori == "v" and dp_pol == "TM":
        a = 1j * norm * q / (k_dp)
        return (-a, -a)

    elif dp_ori == "hx" and dp_pol == "TE":
        a = 1j * norm * np.sin(phi)
        return (-a, -a)

    elif dp_ori == "hx" and dp_pol == "TM":
        a = 1j * norm * k_z / k_dp * np.cos(phi)
        return (a, -a)

    elif dp_ori == "hy" and dp_pol == "TE":
        a = 1j * norm * np.cos(phi)
        return (a, a)

    elif dp_ori == "hy" and dp_pol == "TM":
        a = 1j * norm * k_z / k_dp * np.sin(phi)
        return (a, -a)

    elif dp_ori == "h" and dp_pol == "TE":
        a = norm / np.sqrt(2.)
        return (a, a)

    elif dp_ori == "h" and dp_pol == "TM":
        a = norm / np.sqrt(2.) * k_z / k_dp
        return (a, -a)

    else:
        raise ValueError("Invalid polarization '%s' or orientation '%s'" %
                         (dp_pol, dp_ori))
