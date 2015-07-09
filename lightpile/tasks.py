# -*- coding: utf-8 -*-
"""
Tasks are what the user starts. Tasks are solved by calling engines
and converting data from user-specified quantities to system quantities the
engines expect.
After the engines are done the results are reconverted to user quantitites and
presented according to the requirements of the tasks (can be ASCII-data or
graphs or both).
"""
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from lightpile.edp import EDP
from lightpile.graphs import Graph_y_a
from quantities import AngularRange, AngularArray, SpectralPoint
from quantities import QunttyValue, QunttyArray
from quantities import SpectralConverter, AngularConverter

__author__ = "Richard Pfeifer"


class DipolestudyTask(object):
    """
    Calculate the emission of a single monochromatic source.

    Results:
        1.) ASCII file to path 'f_a_data_path' with figure data
        2.) ASCII file to path 'p_a_data_path' with figure data
        3.) figure 'graph_f_a' f(angular_quantity)
        4.) figure 'graph_p_a' p(angular_quantity)

        the following data is presented:
        f   ang     array   angularrange_user_quantity
            Re(w)   array
            Im(w)   array
            f       array
            f_hTE   array
            f_hTM   array
            f_vTM   array
            P/P0    float

         p  ang
            Re(w)
            Im(w)
            p
            p_hTE
            p_hTM
            p_vTM
            P/P0

    """

    def __init__(self, stack, dipolestudy, **kwargs):
        """
        Needs a 'stack' object and a 'dipolestudy' configuration section.

        A 'graph_f_a' section is optional in the case the default configuration
        for the graph_f_a is to be changed.
        A 'graph_p_a' section is optional in the case the default configuration
        for the graph_p_a is to be changed.
        """
        self.stack = stack
        self.spectralpoint_user = SpectralPoint(**dipolestudy["spectralpoint"])
        self.angularrange_user = AngularRange(**dipolestudy["angularrange"])
        self.text_f_filename = "dipolestudy_data_f.txt"
        self.text_p_filename = "dipolestudy_data_p.txt"
        self.graph_f_filename = "dipolestudy_graph_f.png"
        self.graph_p_filename = "dipolestudy_graph_p.png"
        self.graph_f_a_userkwargs = {}
        self.graph_p_a_userkwargs = {}
        if 'graph_f_a' in kwargs:
            self.graph_f_a_userkwargs.update(kwargs['graph_f_a'])
        if 'graph_p_a' in kwargs:
            self.graph_p_a_userkwargs.update(kwargs['graph_p_a'])
        self.engine = 'edp'     # 'edp' is the only engine at the moment

    def run(self):
        """
        Executes the job 'dipolestudy'.

        The emission of a single monochromatic source is calculated using the
        LOE-engine at a single wavelength.
        Output is available both text-based an as a graph.

        The user given quantities are converted to system used quantities.
        The LOE-engine is initialized, run and the output generated.
        """

        self.sp_cv = SpectralConverter(self.spectralpoint_user)
        self.sp_sys = self.sp_cv.u_to_s(self.spectralpoint_user)

        self.ang_conv = AngularConverter(self.angularrange_user, self.sp_sys)
        self.ar_sys = self.ang_conv.u_to_s(self.angularrange_user)

        if self.engine == 'edp':
            self.run_EDP()
        else:
            raise ValueError()

        self.print_reports()
        self.plot_f_a()
        self.plot_p_a()

    def run_EDP(self):
        edp = EDP(self.stack, self.sp_sys.value, self.ar_sys)
        edp.calc()

        self.ang = AngularArray(self.angularrange_user.quantity,
                                self.angularrange_user.unit,
                                self.ang_conv.num_s_to_u(edp.us))
        ang_qtty = self.angularrange_user.quantity
        jcb_u_to_ang = self.ang_conv.jcb_s_to_u(edp.us)

        self.f_ang = QunttyArray("emission propability density", "?",
                                 "f({0})".format(ang_qtty),
                                 edp.f_u * jcb_u_to_ang)
        self.f_hTE_ang = QunttyArray("emission propability density hTE", "?",
                                     "f_hTE({0})".format(ang_qtty),
                                     edp.f_hTE_u * jcb_u_to_ang)
        self.f_hTM_ang = QunttyArray("emission propability density hTM", "?",
                                     "f_hTM({0})".format(ang_qtty),
                                     edp.f_hTM_u * jcb_u_to_ang)
        self.f_vTM_ang = QunttyArray("emission propability density vTM", "?",
                                     "f_vTM({0})".format(ang_qtty),
                                     edp.f_vTM_u * jcb_u_to_ang)

        self.p_ang = QunttyArray("power density", "?",
                                 "p({0})".format(ang_qtty),
                                 edp.p_u * jcb_u_to_ang)
        self.p_hTE_ang = QunttyArray("power density hTE", "?",
                                     "p_hTE({0})".format(ang_qtty),
                                     edp.p_hTE_u * jcb_u_to_ang)
        self.p_hTM_ang = QunttyArray("power density hTM", "?",
                                     "p_hTM({0})".format(ang_qtty),
                                     edp.p_hTM_u * jcb_u_to_ang)
        self.p_vTM_ang = QunttyArray("power density vTM", "?",
                                     "p_vTM({0})".format(ang_qtty),
                                     edp.p_vTM_u * jcb_u_to_ang)

        self.Re_w = QunttyArray("Re(w)", "None", "Re(w)", edp.we_u.real)
        self.Im_w = QunttyArray("Im(w)", "None", "Im(w)", edp.we_u.imag)
        self.PinP0 = QunttyValue("emitted power", "P0", "P", edp.P)

    def print_reports(self):
        """
        Print data of emitter / emitting layer in ascii-form.
        """
        print_dipolestudy_report(
            self.stack.emitter, self.spectralpoint_user,
            self.angularrange_user, self.ang,
            self.Re_w, self.Im_w, self.f_ang,
            self.f_hTE_ang, self.f_hTM_ang, self.f_vTM_ang,
            self.PinP0, self.text_f_filename)
        print_dipolestudy_report(
            self.stack.emitter, self.spectralpoint_user,
            self.angularrange_user, self.ang,
            self.Re_w, self.Im_w, self.p_ang,
            self.p_hTE_ang, self.p_hTM_ang, self.p_vTM_ang,
            self.PinP0, self.text_p_filename)

    def plot_f_a(self):

        # Default scale of x-axis is the scale of the angularrange given by
        # the user.
        # If the user didn't specifiy an angularrange scale and it is set to
        # 'None', we set the graph_scale to 'linear'.
        if "xscale" not in self.graph_f_a_userkwargs:
            if self.angularrange_user.scale is not None:
                self.graph_f_a_userkwargs["xscale"] = \
                    self.angularrange_user.scale
            else:
                self.graph_f_a_userkwargs["xscale"] = 'linear'
        # Default scale of y-axis is 'linear'
        if "yscale" not in self.graph_f_a_userkwargs:
            self.graph_f_a_userkwargs["yscale"] = 'linear'

        # print(self.graph_f_a_userkwargs)
        g = Graph_y_a(self.ang, self.f_ang, self.f_hTE_ang, self.f_hTM_ang,
                      self.f_vTM_ang, self.spectralpoint_user,
                      self.graph_f_filename,
                      **self.graph_f_a_userkwargs)
        g.plot()

    def plot_p_a(self):
        # Default scale of x-axis is the scale of the angularrange given by
        # the user.
        # If the user didn't specifiy an angularrange scale and it is set to
        # 'None', we set the graph_scale to 'linear'.
        if "xscale" not in self.graph_p_a_userkwargs:
            if self.angularrange_user.scale is not None:
                self.graph_p_a_userkwargs["xscale"] = \
                    self.angularrange_user.scale
            else:
                self.graph_p_a_userkwargs["xscale"] = 'linear'
        # Default scale of y-axis is 'linear'
        if "yscale" not in self.graph_p_a_userkwargs:
            self.graph_p_a_userkwargs["yscale"] = 'linear'

        g = Graph_y_a(self.ang, self.p_ang, self.p_hTE_ang, self.p_hTM_ang,
                      self.p_vTM_ang, self.spectralpoint_user,
                      self.graph_p_filename,
                      **self.graph_p_a_userkwargs)
        g.plot()


def print_dipolestudy_report(emitter, spectralpoint,
                             angularrange_user, angulardata, re_w, im_w,
                             angpwr_quantity_total, angpwr_quantity_hTE,
                             angpwr_quantity_hTM, angpwr_quantity_vTM,
                             power, filename):
    """
    Print data of emitter / emitting layer in ascii-form.
    """

    # header
    header_title = (
        "{0} {1} for dipole emitter\n"
        "-------------------------------------------------\n\n"
        .format(angpwr_quantity_total.quantity, angpwr_quantity_total.symbol))
    orientation = (
        " emitter orientation distribution: horizontal {0:.3f}"
        "  vertical {1:.3f}\n"
        .format(emitter.orientation_weights[0],
                emitter.orientation_weights[1]))
    spectral = (
        " {0:s}: {1:.1f}{2:s}\n".format(
            spectralpoint.quantity, spectralpoint.value, spectralpoint.unit))
    angular = " {0:s}: {1:.3f} to {2:.3f}, {3} points ({4})".format(
        angularrange_user.quantity,
        angularrange_user.begin,
        angularrange_user.end,
        angularrange_user.sampling_points,
        angularrange_user.scale)

    header_parameters = "".join([
        "Parameters\n", orientation, spectral, angular, "\n"])

    header_results = "Results\n {0} {1} = {2:.6f} in {3}\n\n".format(
        power.quantity, power.symbol, power.value, power.unit)

    ang_qty = angularrange_user.quantity

    header_quantities = (
        "{0:14s}{1:16s}{2:16s}{3:16s}{4:16s}{5:16s}{6:16s}").format(
            ang_qty,
            angpwr_quantity_total.symbol,
            re_w.symbol,
            im_w.symbol,
            angpwr_quantity_hTE.symbol,
            angpwr_quantity_hTM.symbol,
            angpwr_quantity_vTM.symbol)

    header = "".join([
        header_title, header_parameters, header_results, header_quantities])

    # angular data
    to_print = np.zeros((angulardata.data.shape[0], 7))

    to_print[:, 0] = angulardata.data
    to_print[:, 1] = angpwr_quantity_total.data
    to_print[:, 2] = re_w.data
    to_print[:, 3] = im_w.data
    to_print[:, 4] = angpwr_quantity_hTE.data
    to_print[:, 5] = angpwr_quantity_hTM.data
    to_print[:, 6] = angpwr_quantity_vTM.data

    np.savetxt(
        filename, to_print, fmt=str('%+.8e'), header=header, comments='# ')
