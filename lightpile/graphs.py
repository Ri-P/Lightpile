# -*- coding: utf-8 -*-
"""
Plot graphs to illustrate results.

Naming scheme:
    _a: A quantity is plotted vs angular data
"""
from __future__ import unicode_literals
from __future__ import print_function

import os
import os.path
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
mpl.rc('font', family='DejaVu Sans')
mpl.rc("text", usetex='false')

__author__ = "Richard Pfeifer"


class Graph_y_a(object):
    """Show y(angulardata) for an emissioncavity at a specific spectral point.

    yscale := [u'linear' | u'log']
    xscale := [u'linear' | u'log']
    xlim_left, xlim_right, ylim_bottom, ylim_top := [None | float]
    plotlineformat := (see matplotlib)
    legendloc := (see matplotlib)
    """

    def __init__(self, angulararray, y_array, y_hTE_array, y_hTM_array,
                 y_vTM_array, spectralpoint, y_ang_filename,
                 xlim_left=None, xlim_right=None,
                 ylim_bottom=None, ylim_top=None,
                 xscale="linear", yscale="log",
                 legendloc='best', lines="total",
                 externaldata=None):

        self.y = y_array
        self.y_hTE = y_hTE_array
        self.y_hTM = y_hTM_array
        self.y_vTM = y_vTM_array
        self.angulararray = angulararray
        self.spectralpoint = spectralpoint

        self.y_ang_filename = y_ang_filename

        def replace_auto(arg):
            if arg == "auto":
                return None
            else:
                return arg

        self.xlim_left = replace_auto(xlim_left)
        self.xlim_right = replace_auto(xlim_right)
        self.ylim_bottom = replace_auto(ylim_bottom)
        self.ylim_top = replace_auto(ylim_top)
        if xscale not in ['linear', 'log']:
            raise ValueError(
                "Error: xscale='{0}' "
                "but needs to be ['linear' | 'log']".format(xscale))
        if yscale not in ['linear', 'log']:
            raise ValueError(
                "Error: yscale='{0}' "
                "but needs to be ['linear' | 'log']".format(yscale))
        self.xscale = xscale
        self.yscale = yscale
        self.legendloc = legendloc
        self.lines = lines
        self.externaldata = externaldata

    def _set_labels(self):
        title = u'{0} for {1:s} = '.format(self.y.symbol,
                                           self.spectralpoint.symbol) + \
                u'{0:.1f}'.format(self.spectralpoint.value) + \
                u'{0:s}'.format(self.spectralpoint.unit)
        self.ax.set_title(title)
        self.ax.set_xlabel(self.angulararray.symbol)
        self.ax.set_ylabel(self.y.symbol)

        self.ax.legend(loc=self.legendloc)

    def _add_externaldata(self):
        """Load y(x)-data from external file and add this as plot to ax.

        Expects to find two column numeric data in a file at absolute path
        self.externaldata["path"] or relative path (relative to lightpile main
        directory).
        """
        if self.externaldata is None:
            return
        ext_path = None
        # assume absolute path
        ext_abspath = os.path.expanduser(
            os.path.normpath(self.externaldata["path"]))
        if os.path.isfile(ext_abspath):
            ext_path = ext_abspath
        # allow path relative to lightpile main directory
        else:
            ext_relpath = os.path.join(os.path.dirname(__file__), os.pardir,
                                       self.externaldata["path"])
            if os.path.isfile(ext_relpath):
                    ext_path = ext_relpath
            else:
                msg = (
                    "Error: File {0} not found at either {1} or {2}.").format(
                        self.externaldata["path"],
                        ext_abspath,
                        ext_relpath)
                raise Exception(msg)
        try:
            with open(ext_path, 'r') as datafile:
                external = np.loadtxt(datafile)
                external_x = external[:, 0]
                external_y = external[:, 1]
        except:
            print("Error loading external data from {0}".format(
                self.externaldata["path"]))
            raise
        self.ax.plot(external_x, external_y, 'k-',
                     label=self.externaldata["label"])

    def plot(self):

        fig = plt.figure()
        self.ax = fig.add_subplot(111)

        if self.lines == "total":
            self.ax.plot(self.angulararray.data, self.y.data,
                         label="{0}".format(self.y.symbol))
        elif self.lines == "all":
            self.ax.plot(self.angulararray.data, self.y.data,
                         label="{0}".format(self.y.symbol))
            self.ax.plot(self.angulararray.data, self.y_hTE.data, 'g--',
                         label="{0}".format(self.y_hTE.symbol))
            self.ax.plot(self.angulararray.data, self.y_hTM.data, 'r--',
                         label="{0}".format(self.y_hTM.symbol))
            self.ax.plot(self.angulararray.data, self.y_vTM.data, 'r-.',
                         label="{0}".format(self.y_vTM.symbol))
        else:
            raise ValueError

        self._add_externaldata()

        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_xlim(left=self.xlim_left, right=self.xlim_right)
        self.ax.set_ylim(bottom=self.ylim_bottom, top=self.ylim_top)

        self._set_labels()

        plt.savefig(self.y_ang_filename)
