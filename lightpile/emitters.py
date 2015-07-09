# -*- coding: utf-8 -*-
"""
Emitters and ensembles of emitters.

The OrientedEmitter is the building block for all kind of ensembles. It is
treated as a single two-level-system with stochastic decay mechanisms and
therefore cannibalistic decay channels (propability for a specific decay
channel depends on the optical environment). Every oriented emitter has an
internal quantum efficiency describing the strength of internal decay for the
unperturbed emitter in vacuum. The real quantum efficiency is determined by
the optical environment.

In contrary the ensembles (SpectralEnsemble, SpatialEnsemble) are collections
of dipole emitters with fixed propabilities (relative Anteilshäufigkeit)
independent of the optical environment.
"""

from __future__ import print_function
from __future__ import unicode_literals

import numpy
import scipy.interpolate.interpolate as interp
from lightpile.material import expand_path, parse_nkfile

__author__ = "Richard Pfeifer"


class EmitterEnsemble(object):

    def __init__(self):
        pass


class SinglePlaneEnsemble(EmitterEnsemble):

    def __init__(self):
        pass


class OrientedEmitter(SinglePlaneEnsemble):
    """
    An emitter with an orientation.

    parameters:
      orientation ('iso' | 'hor' | 'vert' | float)
                  A float in the range [0, 1] is interpreted as *horizontal*
                  weight
      iqe  float   Internal quantum efficiency
    """
    def __init__(self, orientation, iqe=1.0):
        self.iqe = iqe
        if orientation in ['iso', 'isotropic']:
            self.orientation_weights = [2./3., 1./3.]
        elif orientation in ['hor', 'horizontal']:
            self.orientation_weights = [1., 0.]
        elif orientation in ['vert', 'vertical']:
            self.orientation_weights = [0., 1.]
        elif ((0 <= orientation) and (orientation <= 1.)):
            self.orientation_weights = [orientation, 1.-orientation]
        else:
            raise ValueError("orientation is invalid.")


class SpectralEnsemble(SinglePlaneEnsemble):
    """
    Description of ensemble of dipoles with spectral distribution.
    """
    def __init__(self, pl_path, oriented_emitter):
        """
        """
        self.wl_o = []
        self.pl_o = []

        pl_path = expand_path(pl_path)
        self.wl_o, self.pl_o = parse_nkfile(pl_path)

        # expect self.pl_o to contain complex data with imaginary part == 0
        if not (numpy.array(self.pl_o).imag == 0).all():
            raise Exception("Expected exactly two data columns in PL-file " +
                            pl_path + " and found three.")
        self.interp_pl = interp.interp1d(
            self.wl_o, self.pl_o, kind='linear',
            copy=False, bounds_error=True)

    def get_PL(self, wavelengths):
        return self.interp_pl(wavelengths)
