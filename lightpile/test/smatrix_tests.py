# -*- coding: utf-8 -*-

import numpy as np
import unittest
import numpy.testing

import lightpile.smatrix as smatrix
import lightpile.calc as calc


class TestConfExtinctiveAir(object):
    """
    A complete test stack and incidence parameters to check s-matrix funcs.

    Consists of six layers of 'air' with absorption as a parameter. So only
    the phase propagation and absorption is testet.
    Incidence is perpendicular (theta=0) so TE and TM cases are equal.
    """

    def __init__(self, k_extinction):
        self.k_extinction = k_extinction
        self.Ns = (complex(1, -k_extinction) *
                   np.ones((6,), dtype="complex64"))
        self.mus = np.ones((6,), dtype="complex64")
        self.hs = np.array([0., 0.25, 0.5, 0.75, 1.0, 0.],
                           dtype="float32")
        self.k_v = 2. * np.pi / 1
        self.k_t = 0.
        self.qs = calc.kz(self.k_v, self.Ns, self.k_t)

        self._calc_expected_stildes()
        self._calc_expected_Ss()
        self._calc_expected_fieldcoefficients()

    def _calc_expected_stildes(self):
        """
        Calculate the expected layer s-matrices for all layers.
        """
        identity = np.identity(2, dtype="complex64")
        self.expected_stildes = [
            identity,
            (complex(0, -1) * np.exp(-self.k_extinction * np.pi / 2.) *
             identity),
            complex(-1, 0) * np.exp(-self.k_extinction * np.pi) * identity,
            (complex(0, 1) * np.exp(-self.k_extinction * 3 * np.pi / 2.) *
             identity),
            np.exp(-self.k_extinction * 2. * np.pi) * identity]

    def _calc_expected_Ss(self):
        """
        Calculate expected stack s-matrices S(n, m).

        With L = 4 there are interfaces 0, 1, 2, 3, 4

        S(i, j) are S(0, 0), S(0, 1), S(0, 2), S(0, 3), S(0, 4)
                             S(1, 1), S(1, 2), S(1, 3), S(1, 4)
                                      S(2, 2), S(2, 3), S(2, 4)
                                               S(3, 3), S(3, 4)
                                                        S(4, 4)
        They are saved as member self.expected_Ss
        which is a dictionary with the tuples (n, m) as keys and the matrices
        S(n, m) as values.
        """
        S = {(0, 0): [complex(1, 0), 0.],
             (0, 1): [complex(0, -1), 0.25],
             (0, 2): [complex(0, 1), 0.75],
             (0, 3): [complex(-1, 0), 1.5],
             (0, 4): [complex(-1, 0), 2.5],
             (1, 1): [complex(0, -1), 0.25],
             (1, 2): [complex(0, 1), 0.75],
             (1, 3): [complex(-1, 0), 1.5],
             (1, 4): [complex(-1, 0), 2.5],
             (2, 2): [complex(-1, 0), 0.5],
             (2, 3): [complex(0, -1), 1.25],
             (2, 4): [complex(0, -1), 2.25],
             (3, 3): [complex(0, 1), 0.75],
             (3, 4): [complex(0, 1), 1.75],
             (4, 4): [complex(1, 0), 1.]}

        # contruct S-Matrices from phase and extinction length
        identity = np.identity(2, dtype="complex64")
        for nm_tuple, phase_h in S.iteritems():
            extinction = np.exp(-self.k_extinction * 2. * np.pi * phase_h[1])
            S[nm_tuple] = phase_h[0] * identity * extinction

        self.expected_Ss = S

    def _calc_expected_fieldcoefficients(self):
        """
        Calculate the expected field coefficients u(i), d(i).

        We are interested in the field coefficients of all layers
        i = [0, ..., L+1] with given coefficients for the incoming fields
        from the halfspaces u(0) and d(L+1).
        Here we chose u(0) = 1. and d(5) = 1.
        """
        l_up = np.array([0., 0., 0.25, 0.75, 1.5, 2.5])
        l_down = np.array([2.5, 2.5, 2.25, 1.75, 1.0, 0.])

        phase_up = np.exp(-1j * 2. * np.pi * l_up)
        phase_down = np.exp(-1j * 2. * np.pi * l_down)
        ext_up = np.exp(-self.k_extinction * 2. * np.pi * l_up)
        ext_down = np.exp(-self.k_extinction * 2. * np.pi * l_down)

        self.expected_us = phase_up * ext_up
        self.expected_ds = phase_down * ext_down

# -----------------------------------------
# Perform tests
# -----------------------------------------


def test_redhefferstar_product():
    """
    Test the redheffer star product C = redheffer(A, B) with the case

    a = | 1  2 |   b = | 5  6 |
        | 3  4 |       | 7  8 |
    and expected result
    c = | -0.384615385 -0.153846154 |
        | +0.846153846 -2.461538462 |
    """
    exp_c = np.array([
        [-0.384615385, -0.153846154],
        [+0.846153846, -2.461538462]])
    a = np.array([[1, 2], [3, 4]])
    b = np.array([[5, 6], [7, 8]])
    calc_c = smatrix._redheffer_star(a, b)
    numpy.testing.assert_allclose(calc_c, exp_c, rtol=3e-07, atol=1e-16)


class Case(unittest.TestCase):

    def test_get_layer_s_matrix_phase(self):
        """
        Test phase propagation through vacuum.

          0   1   2    3        4     5   layer-index
            |  |    |      |        |
            |  |    |      |        |
             thickness:  1 quarter wavelength
                         2 half wavelength
                         3 three quarter wavelength
                         4 one wavelength
        With N = (1.-0*j) for all layers we expect s_tilde(p) to be
            s_tilde(p) = | phi  0  |
                         | 0   phi |
        with phi = exp(-j q h).
        """
        k_v = 2. * np.pi / 1.
        k_t = 0.
        Ns = np.ones((6,), dtype="complex64")
        qs = np.sqrt((k_v * Ns)**2 - k_t**2)
        qs = np.where(qs.imag <= 0, qs, -qs)        # here: qs = k_v
        mus = np.ones((6,), dtype="complex64")
        hs = np.array([0., 0.25, 0.5, 0.75, 1., 0.], dtype="float32")
        identity = np.identity(2, dtype="complex64")
        expected_stildes = [identity,
                            complex(0, -1) * identity,
                            complex(-1, 0) * identity,
                            complex(0, 1) * identity,
                            identity]
        for pol in ['TE', 'TM']:
            for p in range(5):
                s_tilde = smatrix.get_layer_s_matrix(p, qs, Ns, mus, hs,
                                                     k_v, pol)
                numpy.testing.assert_allclose(s_tilde, expected_stildes[p],
                                              rtol=3e-07, atol=1e-16)

    def test_get_layer_s_matrix_absorption(self):
        """
        Test s_tilde(p) for absorption in stack of equal materials.

        rtol=3e-7 apparantly is needed because we use float32 and complex64.
        """
        conf_1 = TestConfExtinctiveAir(k_extinction=0.01)

        for pol in ['TE', 'TM']:
            for p in range(5):
                s_tilde = smatrix.get_layer_s_matrix(
                    p, conf_1.qs, conf_1.Ns, conf_1.mus,
                    conf_1.hs, conf_1.k_v, pol)
                numpy.testing.assert_allclose(s_tilde,
                                              conf_1.expected_stildes[p],
                                              rtol=3e-07, atol=1e-16)

    def test_get_layer_s_matrix_exp_comparison(self):
        """
        Compare the explicit function with the implicit function to
        calculate the layer-s-matrix.
        """

        conf_1 = TestConfExtinctiveAir(k_extinction=0.01)

        for pol in ['TE', 'TM']:
            for p in range(len(conf_1.hs) - 1):
                s_exp = smatrix.get_layer_s_matrix_exp(
                    p, conf_1.qs, conf_1.Ns, conf_1.mus, conf_1.hs, pol)
                s_imp = smatrix.get_layer_s_matrix_imp(
                    p, conf_1.qs, conf_1.Ns, conf_1.mus, conf_1.hs,
                    conf_1.k_v, pol)
                numpy.testing.assert_allclose(s_exp, s_imp,
                                              rtol=3e-07, atol=1e-16)

    def test_get_stack_s_matrix(self):
        """
        Test for correct stack S-Matrix with extincitive air-stack.
        """

        conf_1 = TestConfExtinctiveAir(k_extinction=0.01)

        for pol in ['TE', 'TM']:
            for n in range(len(conf_1.hs) - 1):
                for m in range(n, len(conf_1.hs) - 1):
                    S = smatrix.get_stack_s_matrix(
                        n, m, conf_1.qs, conf_1.Ns, conf_1.mus,
                        conf_1.hs, conf_1.k_v, pol)
                    S_expected = conf_1.expected_Ss[(n, m)]
                    numpy.testing.assert_allclose(S_expected, S,
                                                  rtol=7e-07, atol=1e-16)

    def test_get_stack_and_split_s_matrices(self):
        """
        Test for correct S(n, m) and Sts(n, m) with extinctive air stack.

        Sts(n, m) is a list with the S-matrices of the splits of the stack
        (for n=0, m=4).
        The list is expected to hold the following lists of tuples:
                [(St(0, 0), St(1, 4)),
                 (St(0, 1), St(2, 4)),
                 (St(0, 2), St(3, 4)),
                 (St(0, 3), St(4, 4))]
        or in generell:
          Sts is expected to hold the (m-n) tuples
               (S(n, p), S(p+1, m)) with p in [n, ..., m-1], for m>n  or
                is empty [] for m==n
        """

        conf_1 = TestConfExtinctiveAir(k_extinction=0.01)
        L = len(conf_1.hs) - 2
        for pol in ['TE', 'TM']:
            for n in range(L + 1):
                for m in range(n, L + 1):
                    print n, m
                    S, Sps = smatrix.get_stack_and_split_s_matrices(
                        n, m, conf_1.qs, conf_1.Ns, conf_1.mus,
                        conf_1.hs, conf_1.k_v, pol)
                    S_expected = conf_1.expected_Ss[(n, m)]
                    numpy.testing.assert_allclose(S, S_expected,
                                                  rtol=7e-07, atol=1e-16)
                    # Pick the expected split S-matrices
                    Sps_expected = [(conf_1.expected_Ss[(n, p)],
                                     conf_1.expected_Ss[(p + 1, m)])
                                    for p in range(n, m)]
                    # Make sure Sps has not more elements than expected
                    numpy.testing.assert_equal(len(Sps),
                                               len(Sps_expected))
                    for p, [Sp_exp_left, Sp_exp_right] in \
                            enumerate(Sps_expected):
                        numpy.testing.assert_allclose(
                            Sp_exp_left, Sps[p][0], rtol=7e-07, atol=1e-16)
                        numpy.testing.assert_allclose(
                            Sp_exp_right, Sps[p][1], rtol=7e-07, atol=1e-16)

    def test_get_field_coefficients(self):
        """Test for correct calculated field coefficients u[i], d[i] in all
        layers i = [0, ..., L+1].
        """
        conf_1 = TestConfExtinctiveAir(k_extinction=0.01)
        for pol in ['TE', 'TM']:
            us, ds = smatrix.get_field_coefficients(
                conf_1.qs, conf_1.Ns, conf_1.mus,
                conf_1.hs, conf_1.k_v, pol, 1., 1.)
            numpy.testing.assert_allclose(
                conf_1.expected_us, us, rtol=7e-07, atol=1e-16)
            numpy.testing.assert_allclose(
                conf_1.expected_ds, ds, rtol=7e-07, atol=1e-16)
