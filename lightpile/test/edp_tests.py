# -*- coding: utf-8 -*-
"""
Test EDP
"""
from __future__ import unicode_literals
from __future__ import print_function
import numpy as np
from nose.tools import assert_raises

from lightpile.edp import EDP
from lightpile.stacks import PlanarWG, SingleEmissionplaneStack, PlanarStack
from lightpile.quantities import AngularRange
from lightpile.material import Material
from lightpile.emitters import OrientedEmitter

__author__ = "Richard Pfeifer"


def stack_ford():
    """
    Get the stack used for Ford1984 Fig. 4.
    """
    air = Material(1.0)
    ag = Material(complex(0.0767+4.366j))
    air_top = PlanarWG(air, 0.)
    air_bot = PlanarWG(air, 0.05)
    ag_bot = PlanarWG(ag, 200)
    emitter = OrientedEmitter(orientation="vert", iqe=1.0)
    stack_ford = SingleEmissionplaneStack([air_top, emitter, air_bot, ag_bot])
    return stack_ford


def test_ford1984_fixed_sampling():
    """
    Compare EDP-calculation to published data of Ford1984 using fixed
    samplingpoints.
    """
    u_begin = 0.1
    u_end = 1e4
    u_sampling = 30
    a_range_ford = AngularRange(
        quantity="u",
        unit="None",
        begin=u_begin,
        end=u_end,
        sampling_points=u_sampling,
        scale="log")
    wl_ford = 633.

    edp = EDP(stack_ford(), wl_ford, a_range_ford)
    edp.calc()

    expected_us = np.logspace(
        np.log10(u_begin), np.log10(u_end), u_sampling)

    # Depending on the used integration rule for presampled data the expeced
    # results change. The choices for integration rules are
    #  'simpson rule' or 'trapezoidal rule'. We use simpson rule.
    expected_fu_simp = np.array([
        2.26251490e-10,  7.48538564e-10, 2.49354914e-09,  8.44135783e-09,
        2.97679978e-08,  1.19471423e-07,  4.12958094e-08,  4.59544802e-09,
        5.02343278e-09,  7.98190669e-09,  1.51456092e-08,  3.13764552e-08,
        6.74498466e-08,  1.46902522e-07,  3.20916627e-07,  6.99828263e-07,
        1.51834593e-06,  3.26581994e-06,  6.93147240e-06,  1.44198875e-05,
        2.91154047e-05,  5.62299028e-05,  1.01641970e-04,  1.66506172e-04,
        2.35615633e-04,  2.68167554e-04,  2.20773058e-04,  1.12272195e-04,
        2.78883745e-05,  2.38642019e-06
    ])

    # expected_fu_trapz = np.array([
    #     2.2209952640e-10,  7.3480208673e-10,  2.4477898665e-09,
    #     8.2864494495e-09,  2.9221722578e-08,  1.1727899292e-07,
    #     4.0537987838e-08,  4.5111162433e-09,  4.9312469264e-09,
    #     7.8354303257e-09,  1.4867671080e-08,  3.0800661873e-08,
    #     6.6212065183e-08,  1.4420670611e-07,  3.1502742945e-07,
    #     6.8698569979e-07,  1.4904825231e-06,  3.2058884131e-06,
    #     6.8042720600e-06,  1.4155266001e-05,  2.8581105859e-05,
    #     5.5198022429e-05,  9.9776726635e-05,  1.6345060860e-04,
    #     2.3129182464e-04,  2.6324636722e-04,  2.1672163426e-04,
    #     1.1021187688e-04,  2.7376590798e-05,  2.3426268963e-06
    # ])

    np.testing.assert_allclose(expected_us, edp.us, rtol=3e-07, atol=1e-16)
    np.testing.assert_allclose(
        expected_fu_simp, edp.f_u, rtol=3e-07, atol=1e-16)


def test_fail_on_non_u_angular_range():
    """
    Fail with ValueError if quantity of angular_range is not 'u'.
    """
    non_u_angular_range = AngularRange(
        quantity='q', unit='1/nm', begin=0, end=100, sampling_points=10,
        scale="linear")
    assert_raises(ValueError, EDP, stack_ford(), 633, non_u_angular_range)


def test_fail_on_non_SingleEmissionplaneStack():
    """
    Fail with TypeError if stack is not SingleEmissionplaneStack.
    """
    air = Material(1.0)
    ag = Material(complex(0.0767+4.366j))
    air_top = PlanarWG(air, 0.)
    air_bot = PlanarWG(air, 0.05)
    ag_bot = PlanarWG(ag, 200)
    planar_stack = PlanarStack([air_top, air_bot, ag_bot])
    angular_range = AngularRange(
        quantity="u", unit="None", begin=0.1, end=1e4, sampling_points=20,
        scale="log"
    )

    assert_raises(TypeError, EDP, planar_stack, 633, angular_range)
