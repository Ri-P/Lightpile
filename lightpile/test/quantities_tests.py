# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import print_function

import numpy as np

from lightpile.quantities import SpectralPoint, SpectralConverter
from lightpile.quantities import AngularRange, AngularArray, AngularConverter

from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_allclose

__author__ = "Richard Pfeifer"


def test_spectralpoint_conversion():
    """
    Test for correct conversion of spectral point u->s and back s->u.
    """
    sp_dict = {"quantity": "wavelength", "unit": "µm", "value": 0.525}
    spec_pnt = SpectralPoint(**sp_dict)

    sp_conv = SpectralConverter(spec_pnt)
    spec_pnt_sys = sp_conv.u_to_s(spec_pnt)
    assert_equal(spec_pnt_sys.quantity, "wavelength")
    assert_equal(spec_pnt_sys.unit, "nm")
    assert_equal(spec_pnt_sys.value, 525)

    spec_pnt_user = sp_conv.s_to_u(spec_pnt_sys)
    assert_equal(spec_pnt_user.quantity, "wavelength")
    assert_equal(spec_pnt_user.unit, "µm")
    assert_equal(spec_pnt_user.value, 0.525)


def test_invalid_angular_range_conversion_with_nosys_spectralpoint():
    """
    Test for ValueError AngularConversion when called with nosys SpectralPoint.
    """
    spec_pnt_nosys = SpectralPoint(
        quantity="wavelength", unit="µm", value="0.525")
    ang_range_user = AngularRange(
        quantity="q", unit="1/nm", begin=0., end=0.012,
        sampling_points=100, scale="linear")
    assert_raises(ValueError, AngularConverter, ang_range_user, spec_pnt_nosys)


def test_angular_range_conversion():
    """
    Test AngularConversion of AngularRange with regular values.

    Conversion of [q in 1/nm] to [u in 'None']
    """
    spec_pnt_sys = SpectralPoint(
        quantity="wavelength", unit="nm", value=525.)
    ang_range_user = AngularRange(
        quantity="q", unit="1/nm", begin=0., end=0.012,
        sampling_points=100, scale="linear")
    ang_conv = AngularConverter(ang_range_user, spec_pnt_sys)
    ang_range_sys = ang_conv.u_to_s(ang_range_user)
    assert_equal(ang_range_sys.begin, 0.)
    assert_allclose(ang_range_sys.end, 1.00267614, rtol=2e-07, atol=1e-16)
    assert_equal(ang_range_sys.quantity, "u")
    assert_equal(ang_range_sys.unit, "None")
    assert_equal(ang_range_sys.sampling_points, 100)
    assert_equal(ang_range_sys.scale, "linear")

    # reverse the conversion
    ang_range_us2 = ang_conv.s_to_u(ang_range_sys)
    assert_allclose(
        ang_range_user.begin, ang_range_us2.begin, rtol=2e-07, atol=1e-16)
    assert_allclose(
        ang_range_user.end, ang_range_us2.end, rtol=2e-07, atol=1e-16)
    assert_equal(ang_range_us2.quantity, "q")
    assert_equal(ang_range_us2.unit, "1/nm")


def test_angular_array_conversion():
    """
    Test AngularConversion of AngularArray with regular values.

    Conversion of [q in 1/nm] to [u in 'None'].
    """
    spec_pnt_sys = SpectralPoint(
        quantity="wavelength", unit="nm", value=525.)
    ang_array_user = AngularArray(
        quantity="q", unit="1/nm", dataarray=np.linspace(0., 0.024, 11))
    ang_conv = AngularConverter(ang_array_user, spec_pnt_sys)
    ang_array_sys = ang_conv.u_to_s(ang_array_user)
    aa_sys_expected = np.array([0., 0.20053523, 0.40107046, 0.60160568,
                                0.80214091, 1.00267614, 1.20321137,
                                1.4037466, 1.60428183, 1.80481705,
                                2.00535228])
    assert_allclose(
        ang_array_sys.data, aa_sys_expected, rtol=2e-07, atol=1e-16)
    assert_equal(ang_array_sys.quantity, "u")
    assert_equal(ang_array_sys.unit, "None")

    # reverse the conversion
    ang_array_us2 = ang_conv.s_to_u(ang_array_sys)
    assert_allclose(
        ang_array_us2.data, ang_array_user.data, rtol=2e-07, atol=1e-16)
    assert_equal(ang_array_us2.quantity, "q")
    assert_equal(ang_array_us2.unit, "1/nm")


def test_angular_range_begin_nonnegative():
    assert_raises(
        ValueError, AngularRange, quantity="q", unit="1/nm", begin=-1.,
        end=0.012, sampling_points=100, scale="linear")


def test_angular_range_end_greater_then_begin():
    assert_raises(
        ValueError, AngularRange, quantity="q", unit="1/nm", begin=0.,
        end=-1., sampling_points=100, scale="linear")


def test_angular_range_log_begin_positive():
    assert_raises(
        ValueError, AngularRange, quantity="q", unit="1/nm", begin=0.,
        end=1., sampling_points=100, scale="log")


def test_print_spectral_point():
    spec_pnt_sys = SpectralPoint(
        quantity="wavelength", unit="nm", value=525.)
    assert_equal("{}".format(spec_pnt_sys), "wavelength λ=525.000 nm")


def test_print_angular_range():
    ang_range = AngularRange(
        quantity="q", unit="1/nm", begin=0., end=0.012,
        sampling_points=100, scale="linear")
    assert_equal(
        "{}".format(ang_range),
        "q q=0.000 .. 0.012 1/nm with 100 points (linear)")


def test_angular_range_linear_getdata():
    ang_range = AngularRange(
        quantity="q", unit="1/nm", begin=0., end=0.012,
        sampling_points=100, scale="linear")
    assert_allclose(
        ang_range.get_array().data, np.linspace(0., 0.012, 100))


def test_angular_array_print():
    ang_array = AngularArray(
        quantity="u", unit="None", dataarray=np.linspace(0, 1., 3))
    assert_equal("{}".format(ang_array), "u u=[ 0.   0.5  1. ] ")
