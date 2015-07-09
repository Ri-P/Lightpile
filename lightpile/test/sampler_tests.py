# -*- coding: utf-8 -*-
"""
Test sampler
"""
from __future__ import unicode_literals
from __future__ import print_function
import numpy as np
import numpy.testing
from nose.tools import assert_equal


from lightpile.sampler import MVFSampler, MVFFixedSampler

__author__ = "Richard Pfeifer"

# Sampling points of the 21-Point Kronrod rule in the intervall [-1, 1]
# Taken from
#   Quadpack - A subroutine Package for Automatic Integration
#   by R.Piessens, E.deDoncker-Kapenga, C.W.Überhuber, D.K.Kahaner
#   Table 2.2
kronrod21_positive = [
    0.9956571630258080807355272,
    0.9739065285171717200779640,
    0.9301574913557082260012071,
    0.8650633666889845107320966,
    0.7808177265864168970637175,
    0.6794095682990244062343273,
    0.5627571346686046833390000,
    0.4333953941292471907992659,
    0.2943928627014601981311266,
    0.1488743389816312108848260]
kronrod21 = [0.]
for p in [1, 3, 5, 7, 9, 0, 2, 4, 6, 8]:
    kronrod21 += [-kronrod21_positive[p], kronrod21_positive[p]]
kronrod21 = np.array(kronrod21, dtype=np.float64)

kronrod21_sorted = numpy.array(
    [-x for x in kronrod21_positive] +
    [0] +
    kronrod21_positive[::-1])


def lorentz_approx(x, e, a):
    """
    Lorentz approximation to delta distribution scaled by constant *a*.
    """
    if (e <= 0):
        raise ValueError
    if (a == 0):
        raise ValueError
    l = e / np.pi / (x**2 + e**2)
    return a * l, l


def degree2polynom(x, a, b, c):
    """
    y = a*x**2 + b*x + c
    """
    return a * x**2 + b * x + c, b*c + c


def test_sample_with_one_intervall():
    """
    Test the position and handling of the 21 sampling points.
    """
    e_test = 0.1
    a_test = 2.
    a = -10
    b = 10

    x_exp = (kronrod21 + 1.) / 2. * (b - a) + a
    lorentz_exp, lorentz_normed_exp = lorentz_approx(x_exp, e_test, a_test)

    sampler = MVFSampler(
        lorentz_approx, (e_test, a_test), (np.float64, np.float64),
        epsrel=1e-4, epsabs=1e-8, n_intervalls_max=1, do_sort=False)
    xs, [lorentz, lorentz_normed] = sampler.eval_interval(a, b)

    print("n_intrvalls = {0}".format(sampler.n_intervalls))
    print("n_eval = {0}".format(sampler.n_samples))
    print("integral_value = {0}".format(sampler.integral_value))

    np.testing.assert_allclose(lorentz, lorentz_exp, rtol=3e-10, atol=1e-10)
    np.testing.assert_allclose(xs, x_exp, rtol=3e-10, atol=1e-10)
    np.testing.assert_allclose(
        lorentz_normed, lorentz_normed_exp, rtol=3e-10, atol=1e-10)
    assert_equal(sampler.n_intervalls, 1)
    assert_equal(sampler.n_samples, 21)


def test_sorting_sample():
    """
    Test the sorting of the sample.
    """
    e_test = 0.1
    a_test = 2.
    a = -10
    b = 10

    x_exp = (kronrod21_sorted + 1.) / 2. * (b - a) + a
    lorentz_exp, lorentz_normed_exp = lorentz_approx(x_exp, e_test, a_test)

    sampler = MVFSampler(
        lorentz_approx, (e_test, a_test), (np.float64, np.float64),
        epsrel=1e-4, epsabs=1e-8, n_intervalls_max=1, do_sort=True)
    xs, [lorentz, lorentz_normed] = sampler.eval_interval(a, b)

    np.testing.assert_allclose(xs, x_exp)
    np.testing.assert_allclose(lorentz, lorentz_exp)


def test_integration_of_lorentz_delta():
    """
    Test the integration of the lorentz delta approximation.
    """
    e_test = 0.001
    a_test = 2.
    a = -100
    b = 100

    i_exp = a_test * (np.arctan(b/e_test) - np.arctan(a/e_test)) / np.pi

    sampler = MVFSampler(
        lorentz_approx, (e_test, a_test), (np.float64, np.float64),
        epsrel=1e-6, epsabs=1e-6, do_sort=False)
    xs, [lorentz, lorentz_normed] = sampler.eval_interval(a, b)

    np.testing.assert_allclose(
        sampler.integral_value, i_exp, rtol=1e-6, atol=1e-6)


def test_fixed_sampling():
    """
    Evaluate and integrate an degree 2 polynom using a list of
    fixed sampling points.
    """
    a, b, c = 3, -2, 1
    x0 = -100
    x1 = 100

    i_exp = (a / 3. * (x1**3 - x0**3) +
             b / 2. * (x1**2 - x0**2) +
             c * (x1 - x0))
    sampler = MVFFixedSampler(
        degree2polynom, (a, b, c), (np.float64, np.float64))
    x_array, y_arrays = sampler.eval_points(np.linspace(x0, x1, 1000))

    np.testing.assert_allclose(
        sampler.integral_value, i_exp, rtol=1e-7, atol=0)
    np.testing.assert_allclose(x_array, np.linspace(x0, x1, 1000))
