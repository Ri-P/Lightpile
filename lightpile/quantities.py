# -*- coding: utf-8 -*-
"""
Handle physical quantities consisting of value and unit.
"""
from __future__ import unicode_literals
from __future__ import print_function
import math
import numpy as np

__author__ = "Richard Pfeifer"


class Quantity(object):
    """
    Abstract base class
    """

    def __init__(self, quantity, unit):
        self.quantity = quantity
        self.unit = unit


class Spectral(Quantity):

    _symbols = {"wavelength": "λ",
                "energy": "E"}

    def __init__(self, quantity, unit):
        super(Spectral, self).__init__(quantity, unit)
        self.symbol = self._symbols[quantity]


class SpectralPoint(Spectral):

    def __init__(self, quantity, unit, value):
        super(SpectralPoint, self).__init__(quantity, unit)
        self.value = value

    def __str__(self):
        return "{0:s} {1:s}={2:.3f} {3:s}".format(
            self.quantity, self.symbol, self.value, self.unit)


class Angular(Quantity):

    _symbols = {"u": "u",
                "q": "q"}   # q = k_transversal

    def __init__(self, quantity, unit):
        super(Angular, self).__init__(quantity, unit)
        self.symbol = self._symbols[quantity]


class AngularRange(Angular):
    """
    Range of angular values with at least *begin* and *end*.

    It can be an
    - explicit range with
        *sampling_points* := integer
        *scale* := 'linear' | 'log'
    - implicit range with
        *sampling_points* := 'adaptive'
        *scale* := None

    Can generate a populated array of the sampling points in case of an
    explicit range, otherwise raises 'TypeError'.
    """
    def __init__(self, quantity, unit, begin, end, sampling_points, scale):
        super(AngularRange, self).__init__(quantity, unit)
        if (begin < 0):
            raise ValueError((
                "Error: 'angular_range.min={0}' needs to be nonnegative."
                ).format(begin))
        if (scale == "log" and begin <= 0):
            raise ValueError((
                "Error: 'angular_range.min={0}' needs to be positive "
                "if angular_range.scale is 'log'."
                ).format(begin))
        if (end <= begin):
            raise ValueError((
                "Error: 'angular_range.max={0}' needs to be larger then "
                "'angular_range.min={1}'."
                ).format(end, begin))

        self.begin = begin
        self.end = end
        self.sampling_points = sampling_points
        self.scale = scale

    def __str__(self):
        if self.unit == "None":
            unit = ""
        else:
            unit = self.unit
        return ("{0:s} {1:s}={2:.3f} .. {3:.3f} {4:s} with "
                "{5:d} points ({6:s})").format(
                    self.quantity, self.symbol, self.begin, self.end,
                    unit, self.sampling_points, self.scale)

    def get_array(self):
        """
        Return AngularArray of n='sampling_points' values scaled by 'scale'.
        """

        if not isinstance(self.sampling_points, (int, long)):
            raise TypeError

        if self.scale == "linear":
            array = np.linspace(self.begin, self.end, self.sampling_points)
        elif self.scale == "log":
            array = np.logspace(np.log10(self.begin), np.log10(self.end),
                                self.sampling_points)
        else:
            raise ValueError()
        return AngularArray(self.quantity, self.unit, array)


class AngularArray(Angular):

    def __init__(self, quantity, unit, dataarray):
        super(AngularArray, self).__init__(quantity, unit)
        self.data = dataarray

    def __str__(self):
        if self.unit == "None":
            unit = ""
        else:
            unit = self.unit
        return ("{0:s} {1:s}={2:s} {3:s}".format(
            self.quantity, self.symbol, self.data.__str__(), unit))


class QunttyValue(Quantity):

    def __init__(self, quantity, unit, symbol, value):
        super(QunttyValue, self).__init__(quantity, unit)
        self.symbol = symbol
        self.value = value


class QunttyArray(Quantity):

    def __init__(self, quantity, unit, symbol, dataarray):
        super(QunttyArray, self).__init__(quantity, unit)
        self.symbol = symbol
        self.data = dataarray


class SpectralConverter(object):
    """
    Converts spectral units from user_quantities to system_quantity.

    System quantity: wavelength in nm
    Supports user quantities:
        'wavelength' in ['nm', 'µm']
    """
    def __init__(self, spectralpoint_user):

        self.user_quantity = spectralpoint_user.quantity
        self.user_unit = spectralpoint_user.unit
        self.sys_quantity = "wavelength"
        self.sys_unit = "nm"

        if self.user_quantity == 'wavelength':
            self.to_user_function = lambda x: x
            if self.user_unit == 'nm':
                self.to_sys_factor = 1.
            elif self.user_unit == 'µm':
                self.to_sys_factor = 1000.
            else:
                raise ValueError
            self.to_user_factor = 1. / self.to_sys_factor
            self.to_sys_function = lambda x: x
        else:
            raise ValueError

    def num_u_to_s(self, numeric_value):
        return self.to_sys_function(self.to_sys_factor * numeric_value)

    def num_s_to_u(self, numeric_value_sys):
        return self.to_user_factor * self.to_user_function(numeric_value_sys)

    def u_to_s(self, sp_user):
        """
        Convert spectral data to wavelength in nm
        """
        if not ((sp_user.quantity == self.user_quantity) and
                (sp_user.unit == self.user_unit)):
            raise ValueError(
                "Error: Expected quantity {0:s} instead of {1:s}."
                "or unit {2:s} instead of {3:s}".format(
                    self.user_quantity, sp_user.quantity,
                    self.user_unit, sp_user.unit)
            )
        return SpectralPoint(self.sys_quantity, self.sys_unit,
                             self.num_u_to_s(sp_user.value))

    def s_to_u(self, sp_sys):
        """
        Convert spectral data from wavelength in nm to user_quantity in
        user_unit.
        """
        if not ((sp_sys.quantity == self.sys_quantity) and
                (sp_sys.unit == self.sys_unit)):
            raise ValueError(
                "Error: Expected quantity {0:s} instead of {1:s} or"
                "unit {2:s} instead of {3:s}.".format(
                    self.sys_quantity, sp_sys.quantity,
                    self.sys_unit, sp_sys.unit)
            )
        return SpectralPoint(
            self.user_quantity, self.user_unit,
            self.num_s_to_u(sp_sys.value))


class AngularConverter(object):
    """
    Converts angular data from user_quantities to system_quantity and back.

    System quantity: u in 'None'
    Supported user quantities:
        'u' in ['None']
        'q' in ['1/nm']         q = k_transversal with k=2pi/wavelength

    """

    def __init__(self, angularrange, spectralpoint_sys):

        self.user_quantity = angularrange.quantity
        self.user_unit = angularrange.unit
        self.sys_quantity = "u"
        self.sys_unit = "None"

        if not (spectralpoint_sys.quantity == "wavelength" and
                spectralpoint_sys.unit == "nm"):
            raise ValueError()
        self.spectralpoint = spectralpoint_sys  # wl in nm

        self._set_conversion()

    def _set_conversion(self):
        """
        Calculates conversion factors and functions based on spectralpoint
        and user quantity and unit. This is part of the initialization.
        """
        k_vac = 2.*math.pi / self.spectralpoint.value
        if self.user_quantity == 'u':
            # u = u                 du/du = 1
            self.to_user_function = lambda x: x
            self.to_user_jcb = lambda x: 1.
            if self.user_unit == 'None':
                self.to_sys_factor = 1.
            else:
                raise ValueError
            self.to_user_factor = 1. / self.to_sys_factor
            self.to_sys_function = lambda x: x
            self.to_sys_jcb = lambda x: 1.

        elif self.user_quantity == 'q':
            # q = u * k_vac         du/dq = 1 / k_vac
            self.to_user_function = lambda x: x * k_vac
            self.to_user_jcb = lambda x: 1. / k_vac
            if self.user_unit == '1/nm':
                self.to_sys_factor = 1.
            else:
                raise ValueError
            self.to_user_factor = 1. / self.to_sys_factor
            self.to_sys_function = lambda x: x / k_vac
            self.to_sys_jcb = lambda x: k_vac
        else:
            raise ValueError

    def num_u_to_s(self, num_user):
        """
        Converts numeric values and arrays without checking.

        Converts the given numeric to system quantity and unit. Expects input
        data to be in user quantity and unit as set during initialization.
        """
        return self.to_sys_function(self.to_sys_factor * num_user)

    def num_s_to_u(self, num_sys):
        """
        Converts numeric values and arrays to user quantity without
        checking.

        Converts the given numeric to user quantity and unit as given during
        initialization. Expects input to be in system quantity and unit. This
        is not checked.
        """
        return self.to_user_factor * self.to_user_function(num_sys)

    def jcb_u_to_s(self, num_user):
        """
        Jacobian of conversion from user numeric values to system quantity.

            jcb = d(user) / d(system)

        Example:  F = Int( f(u) du)             u = kt / kv, du/d(kt) = 1/kv
                    = Int( f(u(kt)) * du/d(kt) d(kt))
                    = Int( f(kt) d(kt))
                in u : trapz(f(u), u)
                in kt: trapz(f(kt), kt)   f(kt) = du/d(kt) * f(u) = 1/kv * f(u)

        """
        # This is the Jacobian of the conversion. It is to be multiplied with
        # the integrand.
        return self.to_sys_jcb(self.to_sys_factor * num_user)

    def jcb_s_to_u(self, num_sys):
        """
        Jacobian of conversion from system numeric values to user quantity.
        """
        return self.to_user_factor * self.to_user_jcb(num_sys)

    def u_to_s(self, angular_user):
        """
        Convert user angular data to system angular data ('u' in 'None')

        angular_user may be: - AngularRange object
                             - AngularArray object
        """
        if not ((angular_user.quantity == self.user_quantity) and
                (angular_user.unit == self.user_unit)):
            msg = (
                "Error: Expected quantity {0:s} instead of {1:s}."
                "or unit {2:s} instead of {3:s}").format(
                    self.user_quantity, angular_user.quantity,
                    self.user_unit, angular_user.unit)
            raise ValueError(msg)

        if (hasattr(angular_user, "begin") and
                hasattr(angular_user, "end")):
            begin_sys = self.num_u_to_s(angular_user.begin)
            end_sys = self.num_u_to_s(angular_user.end)

            angular_sys = AngularRange(self.sys_quantity,
                                       self.sys_unit, begin_sys, end_sys,
                                       angular_user.sampling_points,
                                       angular_user.scale)
        elif (hasattr(angular_user, "data")):
            data_sys = self.num_u_to_s(angular_user.data)
            angular_sys = AngularArray(self.sys_quantity,
                                       self.sys_unit, data_sys)
        else:
            raise TypeError()
        return angular_sys

    def s_to_u(self, angular_sys):
        """
        Convert system angular data from system_quantity 'u' in 'None' to
        user_quantity in user_unit.

        angular_user may be: - AngularRange object
                             - AngularArray object
        """

        if not ((angular_sys.quantity == self.sys_quantity) and
                (angular_sys.unit == self.sys_unit)):
            msg = (
                "Error: Expected quantity {0:s} instead of {1:s}."
                "or unit {2:s} instead of {3:s}").format(
                    self.sys_quantity, angular_sys.quantity,
                    self.sys_unit, angular_sys.unit)
            raise ValueError(msg)

        if (hasattr(angular_sys, "begin") and
                hasattr(angular_sys, "end")):
            begin_user = self.num_s_to_u(angular_sys.begin)
            end_user = self.num_s_to_u(angular_sys.end)

            angular_user = AngularRange(self.user_quantity,
                                        self.user_unit, begin_user, end_user,
                                        angular_sys.sampling_points,
                                        angular_sys.scale)
        elif hasattr(angular_sys, "data"):
            data_user = self.num_s_to_u(angular_sys.data)
            angular_user = AngularArray(self.user_quantity, self.user_unit,
                                        data_user)
        else:
            raise TypeError()
        return angular_user
