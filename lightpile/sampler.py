# -*- coding: utf-8 -*-
"""
Adaptivly sample and integrate multivalued functions.
"""
from __future__ import unicode_literals
from __future__ import print_function
import numpy as np
from scipy.integrate import quad, simps

__author__ = "Richard Pfeifer"


def _wfunc(u, func, fargs, cache_r, cache_u):
    """
    Wrapper for *func* to log sampling points and additional results.
    """
    cache_u.append(u)
    r = func(u, *fargs)
    for i_r in xrange(len(cache_r)):
        try:
            cache_r[i_r].append(r[i_r])
        except IndexError:
            print(
                "Number of results of *func* does not match given number of"
                " 'result_types'.")
            raise
    return r[0]


def _create_arrays_from_cache(cache_u, cache_r, result_types):
    """
    Create np.arrays from lists
    """
    result_arrays = []
    for i_r, rtype in enumerate(result_types):
        a = np.array(cache_r[i_r], dtype=rtype)
        result_arrays.append(a)
    sampling_array = np.array(cache_u)
    return sampling_array, result_arrays


def _sort_arrays(sampling_array, result_arrays):
    """
    Rearrange *sampling_array* and all arrays in *result_arrays* using the sort
    sequence of the *sampling_array*.
    """
    sort_index = np.argsort(sampling_array)
    sampling_array = sampling_array[sort_index]
    for i_result in xrange(len(result_arrays)):
        result_arrays[i_result] = result_arrays[i_result][sort_index]
    return sampling_array, result_arrays


class MVFSampler(object):
    """
    Adaptively sample and integrate multivalued functions.

    A function y0, y1, ..., yn = func(x, fargs) is sampled and the points
    x_i, y0_i = func(x_i, fargs), y1_i, ..., yn_i are returned.

    The sampling points x_i are chosen by an adaptive integration routine.
    Therefore *func* is integrated along the first result dimension.

    Uses the scipy.integrate.quadpack algorithm.
    """
    def __init__(self, func, fargs, result_types, epsabs=1.49e-08,
                 epsrel=1.49e-08, n_intervalls_max=50, do_sort=True):
        self.func = func
        self.fargs = fargs
        self.result_types = result_types
        self.epsrel = epsrel
        self.epsabs = epsabs
        self.n_intervalls_max = n_intervalls_max
        self.do_sort = do_sort

    def eval_interval(self, a, b):
        """
        Determines sampling points in [a, b] and returns function values and
        the sampling points.

        Arguments:
          a, b    float   with [a, b] in the domain of *func*

        Returns:
          sampling_array, list_of_result_arrays

          sampling_array ... np.array of sampling points in intervall [a, b]
          list_of_result_arrays ... List of np.arrays of dtype given in
                                    MVFSampler argument *result_types*
        """
        cache_u = []
        cache_r = [[] for rtype in self.result_types]

        quadresult = quad(
            _wfunc, a, b, args=(self.func, self.fargs, cache_r, cache_u),
            full_output=1, limit=self.n_intervalls_max)

        try:
            i, quaderr, info = quadresult
        except ValueError:
            i, quaderr, info, m = quadresult
        self.n_intervalls = info["last"]
        self.n_samples = info["neval"]
        self.integral_value = i

        sampling_array, result_arrays = _create_arrays_from_cache(
            cache_u, cache_r, self.result_types)

        if self.do_sort:
            sampling_array, result_arrays = _sort_arrays(
                sampling_array, result_arrays)

        return sampling_array, result_arrays


class MVFFixedSampler(object):
    """
    Sample and integrate multivalued functions at fixed given points.

    A function y0, y1, ..., yn = func(x, fargs) is sampled and the points
    x_i, y0_i = func(x_i, fargs), y1_i, ..., yn_i are returned.

    The sampling points x_i are provided by the caller.
    """
    def __init__(self, func, fargs, result_types):
        self.func = func
        self.fargs = fargs
        self.result_types = result_types

    def eval_points(self, sorted_sample_points):
        """
        Evaluates *func* at *sorted_sample_points* and returns function values,
        sampling points, and an approximation to the integral of *func*
        along the first result dimension (y0).

        Arguments:
          sorted_sample_points    float   From the domain of *func*

        Returns:
          sampling_array, list_of_result_arrays

          sampling_array ... np.array of sampling points in intervall [a, b]
          list_of_result_arrays ... List of np.arrays of dtype given in
                                    MVFSampler argument *result_types*
        """
        cache_u = []
        cache_r = [[] for rtype in self.result_types]

        for u in sorted_sample_points:
            _wfunc(u, self.func, self.fargs, cache_r, cache_u)

        sampling_array, result_arrays = _create_arrays_from_cache(
            cache_u, cache_r, self.result_types
        )

        self.n_samples = len(cache_u)
        self.integral_value = simps(
            result_arrays[0],
            sampling_array
        )

        return sampling_array, result_arrays
