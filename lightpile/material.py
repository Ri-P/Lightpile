# -*- coding: utf-8 -*-
"""
Build material objects from datafiles.

Convention of refractive index used:
There are two different conventions for the complex refractive index
  a) N = n + i k      with k > 0 in case of losses
  b) N = n - i k      with k > 0 in case of losses
Lightpile uses convention a) with k>=0 to the *user side* and
               convention b) with k>=0 at the *system side*.

For input, both by nk-text-files or direct user input, lightpile expects k>=0
and N = n + i k.
The internally used complex refractive index will then be (as given in
case b)   N = (n, -k)

The conversion of user nk-data to system nk-data happens in the __init__
of materials()
"""
from __future__ import print_function
from __future__ import unicode_literals

import os.path as path
import numpy as np
import scipy.interpolate.interpolate as interp

__author__ = "Richard Pfeifer"


def parse_nkfile(nk_path):
    """
    Reads nk-data from a textfile with nk-data in suitable format.

    Args:
        nk_path (string): valid file path to well-formatted nk-datafile

    Returns:
        wavelengths   list of real values: wavelengths in nm
        ns            list of complex values: (n, +k)
    """

    # prepare reading file
    with open(nk_path, 'r') as f:
        wavelengths = []
        ns = []
        foundThreeColumnFormatLine = False
        linesNotThreeColumns = []

        def construct_wl_n_data(a, b, c):
            wavelengths.append(float(a))
            ns.append(complex(float(b), -float(c)))

        # read file line by line and check format of each line
        for line in f:
            # everything behind '#' is a comment -> ignore it
            line = (line.split('#')[0]).strip()
            if len(line) == 0:
                continue
            try:
                # try if line has three-column format
                a, b, c = line.split()
            except ValueError:
                # line does not have three-column format
                try:
                    a, b = line.split()
                except ValueError:
                    # save line to later count number of these lines and
                    # evaluate them
                    linesNotThreeColumns.append(line)
                else:
                    # we only have two columns
                    # assume k=0
                    c = 0.
                    foundThreeColumnFormatLine = True
                    construct_wl_n_data(a, b, c)
            else:
                # line has three-column format:
                # assume: wavelength in nm      n       k
                foundThreeColumnFormatLine = True
                construct_wl_n_data(a, b, c)

    # check if we understood the data

    # a) we found what we expected (three column format)
    if (foundThreeColumnFormatLine and (len(linesNotThreeColumns) == 0)):
        # asume format like: lambda in nm, n,k
        pass

    # if we did not understand the data, check if it is
    # this special variant
    elif ((not foundThreeColumnFormatLine) and
          (len(linesNotThreeColumns) == 3)):
        # assume format like
        # first line: lambdas in nm separated by white space
        # second line: n separated by white space
        # third line: k separated by white space
        wavelengthInNm = linesNotThreeColumns[0].split()
        wavelengths = map(lambda x: float(x), wavelengthInNm)
        n = linesNotThreeColumns[1].split()
        k = linesNotThreeColumns[2].split()
        ns = map(lambda x, y: complex(float(x), float(y)), n, k)

    # or check this other special variant
    elif ((not foundThreeColumnFormatLine) and
          (len(linesNotThreeColumns) == 2)):
        # assume format like
        # first line: lambdas in nm separated by white space
        # second line: n separated by white space (no absorption)
        wavelengthInNm = linesNotThreeColumns[0].split()
        wavelengths = map(lambda x: float(x), wavelengthInNm)
        n = linesNotThreeColumns[1].split()
        ns = map(lambda x: complex(float(x), 0.), n)

    else:
        # can't figure out format
        print("foundThreeColumnFormatLine {0}".format(
            foundThreeColumnFormatLine))
        print("len(linesNotThreeColumns)= {0}".format(
            len(linesNotThreeColumns)))
        raise Exception("Can't read your nk-file: {0}".format(nk_path))

    return wavelengths, ns


def expand_path(filepath):
    # check file
    try:
        expanded_path = path.expanduser(path.normpath(filepath))
    except AttributeError:
        raise Exception("Is this a valid path: {0}?".format(expanded_path))
    if not path.exists(expanded_path):
        raise Exception("File " + expanded_path + " does not exist.")
    if not path.isfile(expanded_path):
        raise Exception("File " + expanded_path + " is not a file.")
    return expanded_path


class Material(object):
    """
    Optical material with nk-data.

    Can be initialized with single complex value (nondispersive materials) or
    a txt-file with nk-data (dispersive materials).

    All given nk-data have to follow the convention N = n + i k (k>=0) and
    are transformed to the convention N = n - i k  (k>=0).
    """

    def __init__(self, nkdata):
        """
        initializes the material
        *nkdata* can be a single real or complex value or a file in an
        appropriate format (see below)
        """

        self.is_dispersive = None   # is n_c wavelength dependent?
        self.n_c = None             # complex refractive index

        # for saving the data given in the datafile "nkdata"
        self.wl_o = []      # original wavelengthdata (from file)
        self.n_c_o = []     # original complex n data (from file)

        # 1. material without dispersion:
        #    nkdata is a complex numer, valid for all wavelength
        # 2. material with dispersion:
        #    nkdata a path to a textfile listing (lambda n k)
        #    for all relevant wavelengths
        try:
            n_c = complex(nkdata)
        except ValueError:
            self.is_dispersive = True
        else:
            self.is_dispersive = False
            self.n_c = n_c.conjugate()
            if self.n_c.imag > 0:
                msg = ("Positive imaginary part of refractive"
                       " index {0} means gain material!").format(nkdata)
                raise Exception(msg)
            if self.n_c.real < 0:
                raise Exception("Negative real part of refractive "
                                "index is not supported.")

        if self.is_dispersive:
            # check file
            nk_path = expand_path(nkdata)
            #  read file
            wl_o, n_c_o = parse_nkfile(nk_path)
            # convert to N = n - i k convention
            self.n_c_o = np.array(n_c_o).conjugate()
            self.wl_o = np.array(wl_o)
            if (self.n_c_o.imag > 0).any():
                raise ValueError(
                    "Positive imaginary part of refractive"
                    " index {0} means gain material!").format(self.n_c_o)
            if (self.n_c_o.real < 0).any():
                raise Exception("Negative real part of refractive "
                                "index is not supported.")

            self.interp_n = interp.interp1d(
                self.wl_o, self.n_c_o, kind='linear',
                copy=False, bounds_error=True)

    def get_complexN(self, wavelengths):

        if not self.is_dispersive:
            return self.n_c
        else:
            return self.interp_n(wavelengths)
