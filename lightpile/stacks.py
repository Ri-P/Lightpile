# -*- coding: utf-8 -*-
"""
Organize materials of certain thickness into stacks."
"""
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
from lightpile.emitters import SinglePlaneEnsemble

__author__ = "Richard Pfeifer"


class StackItemError(Exception):
    """
    The items in the item_list do not allow the construction of a valid
    stack.
    """
    pass


class Waveguide(object):
    """
    A z-independent optical layer.
    """
    pass


class PlanarWG(Waveguide):
    """
    A x-y-infinite waveguide of isotropic material with a
    fixed thickness.
    """
    def __init__(self, material, thickness):
        self.material = material
        if (thickness < 0):
            raise ValueError(
                "Invalid layer thickness {0:.2f}.".format(thickness))
        self.d = thickness


class Stack(object):
    """
    A list of waveguides. Each waveguide is z-independent.

    L    number of in-between-layers (thin films)
    Lp1  number of interfaces (L+1)
    """

    def __init__(self, wg_list, do_reverse=True):

        wg_list.reverse()
        self.waveguides = wg_list
        self.L = len(wg_list) - 2


class PlanarStack(Stack):
    """
    A list of planar waveguides.
    """

    def __init__(self, item_list):

        wgp_list = self._check_items(item_list)
        super(PlanarStack, self).__init__(wgp_list)

    def _check_items(self, item_list):
        """
        All items have to be instances of PlanarWG. At least two items.
        """

        number_PlanarWG_items = sum(
            [isinstance(item, PlanarWG) for item in item_list])

        if not number_PlanarWG_items >= 2:
            raise StackItemError("Need at least two layers.")

        if number_PlanarWG_items != len(item_list):
            raise StackItemError("All stack items need to be planar layers.")

        return item_list

    def _get_N_array(self, planarStack, wavelengths):
        n_wl = len(wavelengths)
        N_array = np.zeros((len(planarStack.waveguides), n_wl),
                           dtype="complex64")
        for i_wg, wg in enumerate(planarStack.waveguides):
            N_array[i_wg] = wg.material.get_complexN(wavelengths)
        return N_array

    def _get_d_array(self, planarStack):
        return np.array([wg.d for wg in planarStack.waveguides],
                        dtype="float32")

    def get_N_array(self, wavelengths):
        """
        Provide 2D array of comlex N for all layers at all wavelengths
        """
        return self._get_N_array(planarStack=self, wavelengths=wavelengths)

    def get_d_array(self):
        return self._get_d_array(self)


class SingleEmissionplaneStack(Stack):
    """
    Planar stack with single emission plane described by two PlanarStacks
    (bot and top).

    Emitters can be 'OrientedDipole' or 'SpectralEnsemble' but there can be
    only one per stack and in one z-plane.
    """

    def __init__(self, item_list):
        """
        Constructor using an 'item_list' containing 'stack description items'
        provided by lightpile.taskcreator.
        """

        top_wg_list, bot_wg_list, emitter = self._group_items(item_list)

        # add emission layer of thickness 0 with material of top-stack
        # (emitter is considered to be at interface but within top-stack)
        # TODO: Make sure this material has k=0 (is non-extinctive)
        top_wg_list.append(PlanarWG(top_wg_list[-1].material, 0.))
        bot_wg_list.insert(0, PlanarWG(top_wg_list[-1].material, 0.))
        self.top = PlanarStack(top_wg_list)
        self.bot = PlanarStack(bot_wg_list)
        self.emitter = emitter

    def _group_items(self, item_list):
        """
        Check item_list for containing
        - at least three elements
        - exactly one SinglePlaneEnsemble, not at the beginning/end of the list
        - rest of items are PlanarWG items
        """
        if len(item_list) < 3:
            raise StackItemError

        # check for position of single emitter
        is_single_plane_emitter_list = [
            isinstance(item, SinglePlaneEnsemble) for item in item_list]
        if sum(is_single_plane_emitter_list) != 1:
            raise StackItemError
        em_index = is_single_plane_emitter_list.index(True)
        if em_index == 0 or em_index == len(item_list) - 2:
            raise StackItemError

        number_PlanarWG_items = sum(
            [isinstance(item, PlanarWG) for item in item_list])
        if not (number_PlanarWG_items == (len(item_list) - 1)):
            raise StackItemError

        return [item_list[0:em_index], item_list[em_index+1:],
                item_list[em_index]]

    def get_top_N_array(self, wavelengths):
        return self._get_N_array(self.top, wavelengths)

    def get_bot_N_array(self, wavelengths):
        return self._get_N_array(self.bot, wavelengths)

    def get_top_d_array(self):
        return self._get_d_array(self.top)

    def get_bot_d_array(self):
        return self._get_d_array(self.top)

    def print_layers(self):
        for wg in self.top.waveguides.__reversed__():
            print("{0}\t{1:.1f}nm".format(wg.material, wg.d))
        print("{0}".format(self.emitter))
        for wg in self.bot.waveguides.__reversed__():
            print("{0}\t{1:.1f}nm".format(wg.material, wg.d))
