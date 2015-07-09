# -*- coding: utf-8 -*-
"""
Calculate the emitted power of a single oriented dipole inside a planar stack
at a specific wavelength. Separate TE and TM contributions are shown and
the internal quantum efficiency in vacuum is taken into consideration.
"""
from __future__ import unicode_literals
from __future__ import print_function

import numpy as np
from numpy.linalg.linalg import LinAlgError

from lightpile.stacks import SingleEmissionplaneStack
from lightpile.sampler import MVFSampler, MVFFixedSampler
import lightpile.smatrix
import lightpile.dipole
import lightpile.calc

__author__ = "Richard Pfeifer"


class EDP(object):
    """
    Calculate emitted dipole power.

    No z-information, no outcoupled power. Single wavelength
    """

    def __init__(self, stack, wl, angular_range):
        self.wl = wl
        self.angular_range = angular_range
        if not ((self.angular_range.quantity == 'u') and
                (self.angular_range.unit == "None")):
            msg = ("Expected 'u' in unit 'None'. Got '{0}' " +
                   "in unit '{1}'").format(self.angular_range.quantity,
                                           self.angular_range.unit)
            raise ValueError(msg)

        if not isinstance(stack, SingleEmissionplaneStack):
            raise TypeError("Expected instance of class 'PlanarOledStack'.")
        self.stack = stack

        self.radiationborder_cheatshift = 1e-7

        self.N_top = self.stack.top.get_N_array(np.array([self.wl]))[:, 0]
        self.N_bot = self.stack.bot.get_N_array(np.array([self.wl]))[:, 0]
        self.mu_top = np.ones((self.stack.top.L + 2))
        self.mu_bot = np.ones((self.stack.bot.L + 2))
        self.hs_top = np.array([wg.d for wg in self.stack.top.waveguides])
        self.hs_bot = np.array([wg.d for wg in self.stack.bot.waveguides])

        if not (self.N_top[0] == self.N_bot[-1]).all():
            # PlanarOledStack instances create artificial 0-thickness
            # emission layers around the emitter (one for the top stack and
            # one for the bottom stack).
            # The material used is the material of the top layer closest
            # to the emitter. Therefore the emitter is always enclosed by one
            # material only.
            # This should also ensure that the emitter is always enclosed by
            # non-extinctive material (TODO)
            raise ValueError()

    def _get_p_func(self):
        """
        Return the function p(u, **pargs)
        """
        def pfunc(u, k_vac, n_e,
                  L_top, L_bot, N_top, N_bot, mu_top, mu_bot, hs_top, hs_bot,
                  ow_h, ow_v, radiationborder_cheatshift):
            # We cheat for the point u = n_e because our algorithm
            # has a singular matrix and some singlarities there.
            # Actually this is a valid point and we should think about
            # the correct formula for it. Later.
            if abs(u - n_e) < radiationborder_cheatshift:
                u = n_e - radiationborder_cheatshift
            w_top = lightpile.calc.w(N_top, u)
            w_bot = lightpile.calc.w(N_bot, u)
            w_e = w_top[0]
            kz_top = w_top * k_vac
            kz_bot = w_bot * k_vac
            get_s = lightpile.smatrix.get_stack_s_matrix
            # r_bot = r_bot_du = S_bot[0,1]
            # r_top = r_top_ud = S_top[1,0]
            # TODO: smatrix-function for only these relevant components
            try:
                rbot_TE = get_s(0, L_bot, kz_bot, N_bot,
                                mu_bot, hs_bot, k_vac, 'TE')[0, 1]
                rbot_TM = get_s(0, L_bot, kz_bot, N_bot,
                                mu_bot, hs_bot, k_vac, 'TM')[0, 1]
                rtop_TE = get_s(0, L_top, kz_top, N_top,
                                mu_top, hs_top, k_vac, 'TE')[1, 0]
                rtop_TM = get_s(0, L_top, kz_top, N_top,
                                mu_top, hs_top, k_vac, 'TM')[1, 0]
            except LinAlgError:
                print(
                    "Error (singular matrix) with\n\tu={0}\n\tw={1}"
                    .format(u, w_e))
                raise
            rho_TE = 1. - rbot_TE * rtop_TE
            rho_TM = 1. - rbot_TM * rtop_TM

            p_hTE = ow_h * 3./4. * (u/w_e *
                                    ((1.+rbot_TE)*(1.+rtop_TE)) / rho_TE).real
            p_hTM = ow_h * 3./4. * (u*w_e/n_e**2 *
                                    ((1.-rbot_TM)*(1.-rtop_TM)) / rho_TM).real
            p_vTM = ow_v * 3./2. * (u**3 / (w_e * n_e**2) *
                                    ((1.+rbot_TM)*(1.+rtop_TM)) / rho_TM).real
            p = p_hTE + p_hTM + p_vTM
            return p, w_e, p_hTE, p_hTM, p_vTM

        return pfunc

    def calc(self):

        k_vac = 2.*np.pi/self.wl
        n_e = self.N_top[0].real

        ow_h = self.stack.emitter.orientation_weights[0]
        ow_v = self.stack.emitter.orientation_weights[1]
        iqe = self.stack.emitter.iqe
        pfunc_args = (
            k_vac, n_e,
            self.stack.top.L, self.stack.bot.L,
            self.N_top, self.N_bot,
            self.mu_top, self.mu_bot,
            self.hs_top, self.hs_bot,
            ow_h, ow_v, self.radiationborder_cheatshift)
        try:
            u_points = (self.angular_range.get_array()).data
        except TypeError:
            # do adaptive sampling
            sampler = MVFSampler(
                self._get_p_func(), pfunc_args,
                (np.float32, np.complex64, np.float32, np.float32, np.float32),
                epsabs=0, epsrel=1e-7, n_intervalls_max=50, do_sort=True
            )
            sampler_result = sampler.eval_interval(
                self.angular_range.begin, self.angular_range.end)
        else:
            # do fixed sampling
            sampler = MVFFixedSampler(
                self._get_p_func(), pfunc_args,
                (np.float32, np.complex64, np.float32, np.float32, np.float32),
            )
            sampler_result = sampler.eval_points(u_points)

        self.us = sampler_result[0]
        self.p_u = sampler_result[1][0]
        self.we_u = sampler_result[1][1]
        self.p_hTE_u = sampler_result[1][2]
        self.p_hTM_u = sampler_result[1][3]
        self.p_vTM_u = sampler_result[1][4]
        self.P = sampler.integral_value

        self.f_u = self.p_u / ((1. - iqe) / iqe + self.P)
        self.f_hTE_u = self.p_hTE_u / ((1. - iqe) / iqe + self.P)
        self.f_hTM_u = self.p_hTM_u / ((1. - iqe) / iqe + self.P)
        self.f_vTM_u = self.p_vTM_u / ((1. - iqe) / iqe + self.P)
