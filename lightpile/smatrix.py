# -*- coding: utf-8 -*-
"""
Use S-Matrix scheme to propagate phase information through planar stacks.

smatrix implements algorithms to propagate phase information through planar
thin film stacks. It uses the S-Matrix scheme over the T-Matrix scheme for
stability reasons.

References:
    reference Li96-1024:
    "Formulation and comparison of two recursive matrix algorithms for
     modeling layered diffraction gratings"
    Lifeng Li
    J. Opt. Soc. Am. A, Vol.13, No.5, May 1996

    reference Kim07
    "Extended scattering-matrix methode for efficient full parallel
     implementation of rigorous coupled-wave analysis"
    Hwi Kim, Il-Min Lee, and Byoungho Lee
    J. Opt. Soc. Am. A, Vol.24, No.8, 2007

Conventions:
    We assume the exp(-iqz) convention for spatial phase propagation and
    therefore exp(+i omega t) temporal propagation. The complex optical
    index is N = n - ik with positive k in case of absorption.

Attributes:
    c (float): Speed of light in vacuum.
    mu_v (float): Vacuum permeability (magnetic constant).
    eps_v (float): Vacuum permittivity (elecric constant).
"""
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
from scipy.linalg import inv

__author__ = "Richard Pfeifer"

# permeability of vacuum
mu_v = 4.*np.pi*1e-7
c = 299792458.
eps_v = 1./(mu_v * c * c)


# k_v   vacuum wave vector = 2 pi / wavelength
# q     z-component of wavevector
# qs    array of many q
# Ns    array of complex refrative indices
# mus   array of complex relative magnetic permeabilities
# hs    array of layer thicknesses (non-cumulative)

def get_layer_s_matrix(p, qs, Ns, mus, hs, k_v, pol):
    """
    Calculate scattering matrix for layer p for given polarization.

    S-Matrix s_tilde(p) is calculated for a given stack with L thin film
    layers, L+1 interfaces and L+2 total layers (including adjacent
    half-spaces).
    s_tilde(p) includes interface effects of interface p and propagation
    effects through layer p. It does not include interface effects of
    interface p+1.

    Args:
        p (int): Layerindex with 0 <= p <= L
        qs (ndarray, complex64, (L+2,)):
            z-component of wavevector of stack layers
        Ns (ndarray, complex64, (L+2,)):
            complex refractive index N of stack layers
        mus (ndarray, complex64, (L+2,)): permeability of stack layers
        hs (ndarray, float32, (L+2,)): height of stack layers in nm
        k_v (float): vacuum wave vector
        pol (str): polarization ('TE'|'TM')
    Returns:
        ndarray of dtype complex and shape (2,2)

    |u(p+1)| = s_tilde(p) |u(p)  |
    |d(p)  |              |d(p+1)|
    """
    return get_layer_s_matrix_exp(p, qs, Ns, mus, hs, pol)


def get_layer_s_matrix_imp(p, qs, Ns, mus, hs, k_v, pol):
    """Implicitely calculate the layer s-matrix.

    Uses the matrix formulation for the interface boundary condition to
    calculate the interface smatrix. This approach is longer than the
    explicit calculation but serves as documentation and a way to crosscheck
    results.
    """
    if p < 0 or p > (len(hs) - 1):  # p (interface) in [0...L]
        raise ValueError
    if pol == 'TE':
        # W-matrix: tangential fields of upward / downward travelling waves
        #           for s-polarization
        #           | Ey_up   Ey_down  |
        #           | Hx_up   Hx_down  |
        #       Ey_up = Ey_down = 1
        #      -Hx_up = Hx_down = q / (mu * k_vac)
        Wp = np.ones((2, 2), dtype=np.complex)      # W matrix of layer 'p'
        Hxp = qs[p] / (mu_v * mus[p] * k_v)
        Wp[1, 0] = -Hxp
        Wp[1, 1] = Hxp

        Wpp = np.ones((2, 2), dtype=np.complex)      # W matrix of layer 'p+1'
        Hxpp = qs[p+1] / (mu_v * mus[p+1] * k_v)
        Wpp[1, 0] = -Hxpp
        Wpp[1, 1] = Hxpp
    elif pol == 'TM':
        # W-matrix: tangential fields of upward / downward travelling waves
        #           for p-polarization
        #           | Ex_up   Ex_down  |
        #           | Hy_up   Hy_down  |
        #       Ex_up = -Ex_down = q / (n k_vac) = cos(theta)
        #       Hy_up =  Hy_down = n / mu

        Wp = np.ones((2, 2), dtype=np.complex)      # W matrix of layer 'p'
        Exp = qs[p] / (k_v * Ns[p])
        Wp[0, 0] = Exp
        Wp[0, 1] = -Exp

        Hyp = Ns[p] / (mu_v * mus[p])
        Wp[1, 0] = Hyp
        Wp[1, 1] = Hyp

        Wpp = np.ones((2, 2), dtype=np.complex)      # W matrix of layer 'p+1'
        Expp = qs[p+1] / (k_v * Ns[p+1])
        Wpp[0, 0] = Expp
        Wpp[0, 1] = -Expp

        Hypp = Ns[p+1] / (mu_v * mus[p+1])
        Wpp[1, 0] = Hypp
        Wpp[1, 1] = Hypp
    else:
        raise ValueError

    t_p = np.dot(inv(Wpp), Wp)                  # see Li96-1024 equation (7)

    s_p = np.ones((2, 2), dtype=np.complex)      # see Li96-1024 equ (14a)
    s_p[0, 0] = t_p[0, 0] - t_p[0, 1] / t_p[1, 1] * t_p[1, 0]
    s_p[0, 1] = t_p[0, 1] / t_p[1, 1]
    s_p[1, 0] = - t_p[1, 0] / t_p[1, 1]
    s_p[1, 1] = 1. / t_p[1, 1]

    # phi_left and phi_right: see Li96-1024 equ (13a)
    # but beware: - Li uses the exp(+iqz) convention
    #             - Li uses lambda+ for q_up and lambda- for q_down = -q_up
    # here we use exp(-iqz) convention and q_up = q, q_down = -q
    phi_left = np.array([[1., 0],
                         [0., np.exp(1j*(-qs[p])*hs[p])]])
    phi_right = np.array([[np.exp(-1j*qs[p]*hs[p]), 0.],
                          [0., 1.]])

    # layer s-matrix s_tilde                    # see Li96-1024 equ (13a)
    st_p = np.dot(phi_left, np.dot(s_p, phi_right))
    return st_p


def get_layer_s_matrix_exp(p, qs, Ns, mus, hs, pol):
    """
    Calculate scattering matrix for layer p using explicit formulas.

    This method is numerically more stable and should be faster than
    'get_layer_s_matrix_imp'. The resulting s-matrix should be identical,
    aside from numerical rounding errors.

    Args:
        p (int): Layerindex with 0 <= p <= L
        qs (ndarray, complex64, (L+2,)):
            z-component of wavevector of stack layers
        Ns (ndarray, complex64, (L+2,)):
            complex refractive index N of stack layers
        mus (ndarray, complex64, (L+2,)): permeability of stack layers
        hs (ndarray, float32, (L+2,)): height of stack layers in nm
        k_v (float): vacuum wave vector
        pol (str): polarization ('TE'|'TM')
    Returns:
        ndarray of dtype complex and shape (2,2)

    |u(p+1)| = s_tilde(p) |u(p)  |
    |d(p)  |              |d(p+1)|
    """
    if p < 0 or p > (len(hs) - 1):  # p (interface) in [0...L]
        raise ValueError
    # Calculate interface Fresnel coefficients for interface s-matrix
    if pol == 'TE':
        denominator = qs[p] + mus[p] / mus[p] * qs[p+1]
        t_uu = 2 * qs[p] / denominator
        r_ud = (qs[p] - mus[p] / mus[p+1] * qs[p]) / denominator
        t_dd = 2 * qs[p+1] * mus[p] / mus[p+1] / denominator
        r_du = -r_ud
    elif pol == 'TM':
        a = qs[p+1] * Ns[p] / Ns[p+1]
        b = mus[p] / mus[p+1] * qs[p] * Ns[p+1] / Ns[p]
        t_uu = 2 * qs[p] / (a + b)
        r_ud = (b - a) / (a + b)
        t_dd = 2 * mus[p] / mus[p+1] * qs[p+1] / (a + b)
        r_du = -r_ud
    else:
        raise ValueError

    # phase propagation term through layer p
    phi = np.exp(-1j * qs[p] * hs[p])

    # use interface fresnel coefficients (interface s-matrix) and phase
    # propagation term to calculate layer s-matrix
    st_p = np.ones((2, 2), dtype="complex64")
    st_p[0, 0] = phi * t_uu
    st_p[1, 0] = phi**2 * r_ud
    st_p[0, 1] = r_du
    st_p[1, 1] = phi * t_dd

    return st_p


def _redheffer_star(a, b):
    """Calculate redheffer star product of the two (2,2)-matrices a and b.

    from   --g-->|    |--k-->|    |--i-->
                 |  a |      |  b |
           <-h---|    |<-l---|    |<-j--

    calc         --g-->|     |--i-->
                       |  c  |
                 <-h---|     |<-j---
    Args:
        a   matrix of shape (2, 2) and dtype="complex64"
        b   matrix of shape (2, 2) and dtype="complex64"
    Returns:
        c   matrix of shape (2, 2) and dtype="complex64" which is the
            redheffer star product of a and b.
    """
    c_00 = b[0, 0] * a[0, 0] / (1. - a[0, 1] * b[1, 0])
    c_01 = b[0, 1] + b[0, 0] * a[0, 1] / (1. - b[1, 0] * a[0, 1]) * b[1, 1]
    c_10 = a[1, 0] + a[1, 1] * b[1, 0] / (1. - a[0, 1] * b[1, 0]) * a[0, 0]
    c_11 = a[1, 1] * b[1, 1] / (1. - b[1, 0] * a[0, 1])
    return np.array([[c_00, c_01],
                     [c_10, c_11]], dtype="complex64")


def get_stack_s_matrix(n, m, qs, Ns, mus, hs, k_v, pol):
    """Calculates scattering matrix S(n, m) for stack of layers (0, L+1) for
    s/p-polarization.

    |u(m+1)| = S(n,m) |u(n)  |
    |d(n)  |          |d(m+1)|

    Args:
        n (int): Index of first interface to cross (0 <= n <= L)
        m (int): Index of last interface to cross (0 <= m <= L, m >= n)
        qs (ndarray, complex64, (L+2,)):
            z-component of wavevector of stack layers
        Ns (ndarray, complex64, (L+2,)):
            complex refractive index N of stack layers
        mus (ndarray, complex64, (L+2,)): permeability of stack layers
        hs (ndarray, float32, (L+2,)): height of stack layers in nm
        k_v (float): vacuum wave vector
        pol (str): polarization ('TE'|'TM')

    Returns:
        S(n, m) ndarray of dtype complex64 and shape (2, 2)
    """

    assert len(qs) == len(mus)
    assert len(qs) == len(hs)
    assert (qs.imag <= 0).all()

    assert n in range(len(qs)-1)
    assert m in range(len(qs)-1)
    assert n <= m

    # list of scattering matrices S(n,p) with p in [n, ..., m]
    # S_start_n_list = []

    S_pm1 = np.eye(2, dtype=np.complex)     # S_p_minus_1
    for p in range(n, m+1):
        st_p = get_layer_s_matrix(p, qs, Ns, mus, hs, k_v, pol)
        S_p = _redheffer_star(S_pm1, st_p)
        S_pm1 = S_p
    return S_p


def get_stack_and_split_s_matrices(n, m, qs, Ns, mus, hs, k_v, pol):
    """Calculate stack S(n, m) matrix and S-Matrices for all splits of stack.

    Args:
        n (int): Index of first interface to cross (0 <= n <= L)
        m (int): Index of last interface to cross (0 <= m <= L, m >= n)
        qs (ndarray, complex64, (L+2,)):
            z-component of wavevector of stack layers
        Ns (ndarray, complex64, (L+2,)):
            complex refractive index N of stack layers
        mus (ndarray, complex64, (L+2,)): permeability of stack layers
        hs (ndarray, float32, (L+2,)): height of stack layers in nm
        k_v (float): vacuum wave vector
        pol (str): polarization ('TE'|'TM')

    Returns:
        S, Sts
        S: S(n, m) ndarray of dtype complex64 and shape (2, 2)
        Ssplits: [(S(n, p), S(p+1, m)), ...] for p in [n, ..., m-1]
                 List of tuples of S-matrices for all possible splits of the
                 stack (n, m). Empty list for m==n.
    """

    assert len(qs) == len(mus)
    assert len(qs) == len(hs)
    assert (qs.imag <= 0).all()

    assert n in range(len(qs)-1)
    assert m in range(len(qs)-1)
    assert n <= m

    Ssplits_left = []
    Ssplits_right = []
    st_ps = []      # stilde[p] (list of layer S-Matrices)

    # start with left side of splits
    S_pm1 = np.eye(2, dtype="complex64")
    for p in range(n, m + 1):
        st_p = get_layer_s_matrix(p, qs, Ns, mus, hs, k_v, pol)
        st_ps.append(st_p)
        S_p = _redheffer_star(S_pm1, st_p)
        S_pm1 = S_p
        if p < m:
            Ssplits_left.append(S_p)
        else:
            S = S_p

    # now do the right side of the splits
    for p in range(n, m):
        # generate S(p + 1, m)
        S_km1 = np.eye(2, dtype="complex64")
        for k in range(p + 1, m + 1):
            # accessing st_ps: List position is not stack position if n>0
            st_k = st_ps[k - n]
            S_k = _redheffer_star(S_km1, st_k)
            S_km1 = S_k
        Ssplits_right.append(S_k)
    return S, zip(Ssplits_left, Ssplits_right)


def get_stack_and_intermediate_s_matrices(n, m, qs, Ns, mus, hs, k_v, pol):
    """Calculates scattering matrices of stack and all splits of the stack.

    |u(m+1)| = S(n,m) |u(n)  |
    |d(n)  |          |d(m+1)|

    Returns:
        S, Sis

        S:   ((2, 2) matrix, dtype=complex64) stack S-Matrix S(n, m)
        Ssi: Dictionary with (i, j) tuples as keys and S(i, j) as values.
             It contains the 2*(m-n-1) split S-matrices
                 S(n, p), S(p+1, m-1) with p in [n, ..., m-2]


    """
    assert len(qs) == len(mus)
    assert len(qs) == len(hs)
    assert (qs.imag <= 0).all()

    assert n in range(len(qs)-1)
    assert m in range(len(qs)-1)
    assert n <= m

    # list of pairs [S_left, S_right] with
    # S_left  = S(n, p-1)
    # S_right = S(p, m)
    Sp_pair_list = []
    st_p_list = []      # list of layer-p-matrices s_tilde(p)
    S_left_list = []
    S_right_list = []
    # S_p_minus_1
    S_pm1 = np.eye(2, dtype=np.complex)
    # get all layer matrices and all S_left scattering matrices
    for p in range(n, m+1):
        st_p = get_layer_s_matrix(p, qs, Ns, mus, hs, k_v, pol)
        st_p_list.append(st_p)
        S_p = np.array([[
            st_p[0, 0] * S_pm1[0, 0] / (1 - S_pm1[0, 1] * st_p[1, 0]),
            st_p[0, 1] + st_p[0, 0] * S_pm1[0, 1] / (
                1 - st_p[1, 0] * S_pm1[0, 1]) * st_p[1, 1]], [
            S_pm1[1, 0] + S_pm1[1, 1] * st_p[1, 0] / (
                1 - S_pm1[0, 1] * st_p[1, 0]) * S_pm1[0, 0],
            S_pm1[1, 1] / (1 - st_p[1, 0] * S_pm1[0, 1]) * st_p[1, 1]]
        ])
        S_left_list.append(S_p)
        S_pm1 = S_p
    # S_left_list = [S(n, n), S(n, n+1), S(n, n+2), ..., S(n, m)]
    # with (m-n)+1 elements
    #
    # get S_right scattering matrices S(p,m) for p in [n+1, m] (p=n -> S(n, m)
    #                                                         not interesting)
    # add S(n, m) as dummy element to simplify indexing
    S_right_list.append(None)
    for p in range(n+1, m+1):
        S_km1 = np.eye(2, dtype=np.complex)
        for k in range(p, m+1):
            st_k = st_p_list[k]
            S_k = np.array([[
                st_k[0, 0] * S_km1[0, 0] / (1 - S_km1[0, 1] * st_k[1, 0]),
                st_k[0, 1] + st_k[0, 0] * S_km1[0, 1] / (
                    1 - st_k[1, 0] * S_km1[0, 1]) * st_k[1, 1]], [
                S_km1[1, 0] + S_km1[1, 1] * st_k[1, 0] / (
                    1 - S_km1[0, 1] * st_k[1, 0]) * S_km1[0, 0],
                S_km1[1, 1] / (1 - st_k[1, 0] * S_km1[0, 1]) * st_k[1, 1]]
            ])
            S_km1 = S_k
        S_right_list.append(S_k)
    # S_right_list = [None, S(n+1, m), S(n+2, m), S(n+3, m), ..., S(m, m)]
    # with (m-n)+1 elements

    Sp_pair_list = []
    Sp_pair_list.append(None)   # dummy element for p=n
    for p in range(n+1, m+1):
        Sp_pair_list.append([S_left_list[p-1], S_right_list[p]])
    S_stack = S_left_list[m-n]
    return S_stack, Sp_pair_list


def get_field_coefficients(qs, Ns, mus, hs, k_v, pol, u_0, d_LpOne):
    """Calculate field coefficients u_p, d_p for all layers p in [1, ..., L].

    u_0     : float u[0]
    d_LpOne : float d[L+1]
    -----
    u_p [d_p] are the amplitudes of the upward [downward] travelling waves
    at the beginning of layer p directly behind the interface:
                amplitude of layer      at z=
      u_1, d_1:         1                 0 + d0             = 0 +
      u_2, d_2:         2                 sum(d0+...+d1)     = z1
      u_L, d_L:         L                 sum(d0+...+d(L-1)) = z(L-1)

         0      z1      z2      z(L-1)   zL     ---> z  position
         |      |       |  ...  |        |
       0 |  1   |   2   |       |   L    |   L+1        layer
      d0 |  d1  |   d2  |       |   dL   |   d(L+1)     thickness


    Returns
    -------
    us, ds : ndarray (float) of length L+2
      a = u | d
        a[1:L+1] = a(p) for layers p in [1, ..., L]


    """

    L = len(qs)-2
    S_stack, Sp_pair_list = get_stack_and_intermediate_s_matrices(
        0, L, qs, Ns, mus, hs, k_v, pol)
    return _field_coefficients(qs, Ns, mus, hs, k_v, pol, u_0, d_LpOne,
                               S_stack, Sp_pair_list)


def _field_coefficients(qs, Ns, mus, hs, k_v, pol, u_0, d_LpOne, S_stack,
                        Sp_pair_list):

    assert (qs.imag <= 0).all()
    L = len(qs)-2
    us = np.zeros((L+2), dtype=np.complex)
    ds = np.zeros((L+2), dtype=np.complex)

    us[0] = u_0
    ds[L+1] = d_LpOne
    us[L+1] = S_stack[0, 0] * u_0 + S_stack[0, 1] * d_LpOne
    ds[0] = S_stack[1, 0] * u_0 + S_stack[1, 1] * d_LpOne

    for p in range(1, L+1):
        # see Kim07 equation (37a), (37b), (37c), (37d) with
        # u_p = C^{(0, L+1)+}_{a,(p)} + C^{(0, L+1)+}_{b, (p)}
        # d_p = C^{(0, L+1)-}_{a,(p)} + C^{(0, L+1)-}_{b, (p)}
        # and
        # | u_p | = | F_11 F_12 | | u_0  |
        # | d_p |   | F_21 F_22 | | d_L+1|
        # with
        # F_11 = (37a)
        # F_12 = (37c)
        # F_21 = (37b)
        # F_22 = (37d)
        #
        # left ->  (0, n-1) -> (b)ottom
        # right -> (n, L+1) -> (t)op

        # S = | t_uu r_ud | = | t_u r_u |
        #     | r_du t_dd |   | r_d t_d |
        # r_bu -> r_bottom_up -> r_ud of S_bottom -> S_left[0,1]
        t_bu = Sp_pair_list[p][0][0, 0]
        r_bu = Sp_pair_list[p][0][0, 1]
        # r_bd = Sp_pair_list[p][0][1, 0]
        # t_bd = Sp_pair_list[p][0][1, 1]

        # t_tu = Sp_pair_list[p][1][0, 0]
        # r_tu = Sp_pair_list[p][1][0, 1]
        r_td = Sp_pair_list[p][1][1, 0]
        t_td = Sp_pair_list[p][1][1, 1]

        F_11 = 1. / (1. - r_bu * r_td) * t_bu
        F_12 = r_bu / (1. - r_td * r_bu) * t_td
        F_21 = r_td / (1. - r_bu * r_td) * t_bu
        F_22 = 1. / (1. - r_td * r_bu) * t_td

        us[p] = F_11 * u_0 + F_12 * d_LpOne
        ds[p] = F_21 * u_0 + F_22 * d_LpOne

    return us, ds


def positions_to_layers(hs, zs):
    """Find the layer of positions zs with hs being the thickness of the layers.

    z is in layer i if and only if
        z_i-1 < z <= z_i
              or
        z in (z_i-1, z_i]
    """
    # test if z < layerboundary for all layers
    # add 'True's to get layerindex of z
    return (zs > hs.cumsum()[:-1, np.newaxis].repeat(len(zs), 1)).sum(axis=0)


def expand_field_coefficients_2D(us, ds, hs, Ns, q, z_expand):
    """See 'expand_field_coefficients' but with 2D arrays 'us' 'ds' with
    first axis being the zlayer-axis.
        us.shape = (n_layer, n_other)
        hs.shape = (n_layer)
        q.shape  = (n_layer)
        z_expand.shape = (n_expand)
    """
    ls = positions_to_layers(hs, z_expand)      # [n_zex]
    q_zex = q[ls].swapaxes(0, 1)                 # [n_other, n_zex]
    N_zex = Ns[ls]
    u0_zex = us[ls].swapaxes(0, 1)               # [n_other, n_zex]
    d0_zex = ds[ls].swapaxes(0, 1)               # [n_other, n_zex]
    z0 = np.cumsum(hs) - hs
    z0_zex = z0[ls]
    z_diff = (z_expand - z0_zex)

    u_phase = np.exp(-1j * q_zex * abs(z_diff))
    d_phase = np.exp(1j * q_zex * (z_diff))

    u_zex = u0_zex * u_phase   # [n_zex, n_other]
    d_zex = d0_zex * d_phase   # [n_zex, n_other]

    return (u_zex.swapaxes(0, 1), d_zex.swapaxes(0, 1),
            q_zex.swapaxes(0, 1), N_zex)


def expand_field_coefficients(us, ds, hs, Ns, qs, z_expand):
    """Calculate the field coefficient at positions z_expand from the
       field coefficients us / ds at the layer interfaces using the propagation
       factors.

       Propagation factors are u_p = e^(-j q z') and d_p = e^(j q z')
       for the given stack positions z'.

       z' = z - z0_i    z0_i is the position of the lower layerborder of
                        layer i (the layer with z in [z0_i, z0_i+1)  )
    """
    # wir wollen aus n_layer Werten n_z Werte berechnen

    assert (qs.imag <= 0).all()

    ls = positions_to_layers(hs, z_expand)
    q_zex = qs[ls]
    N_zex = Ns[ls]
    u0_zex = us[ls]
    d0_zex = ds[ls]
    z0 = np.cumsum(hs)-hs       # z0 = z_p - h_p for layer p
    z0_zex = z0[ls]

    u_zex = u0_zex * np.exp(-1j * q_zex * abs(z_expand - z0_zex))
    d_zex = d0_zex * np.exp(1j * q_zex * (z_expand - z0_zex))

    return u_zex, d_zex, q_zex, N_zex
