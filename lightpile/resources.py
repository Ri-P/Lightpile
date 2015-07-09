# -*- coding: utf-8 -*-
"""
Find resources like nk- or pl-files.
"""
from __future__ import unicode_literals
from __future__ import print_function
import os
import os.path

__author__ = "Richard Pfeifer"


def get_resource_path(pathstring):
    """Use *pathstring* to locate resource and return absolute path.

    Accepts and searches in descending order:
            - absolute path (userdirectory shortener allowed)
            - path relative to current working directory
    TODO: Add configuration to allow user directory which is added to
          search path for resources.
    """

    searchdirs = [os.getcwd()]
    failed_paths = []
    resource_path = None

    # expand userdir shortener
    try:
        expanded_path = os.path.expanduser(pathstring)
    except AttributeError:
        raise Exception("Error, not a valid path: {0}.".format(pathstring))

    if os.path.isabs(expanded_path):
        if os.path.isfile(expanded_path):
            resource_path = expanded_path
        else:
            failed_paths.append(expanded_path)
    else:
        try_path = expanded_path
        for searchdir in searchdirs:
            try_path = os.path.realpath(os.path.join(searchdir, expanded_path))
            if os.path.isfile(try_path):
                resource_path = try_path
                break
            else:
                failed_paths.append(try_path)
    # If unsuccessfull, announce failure and tried paths for user learning
    # process.
    if None == resource_path:
        raise Exception(
            "Error. Found no file {0} although I searched at {1}."
            .format(pathstring, failed_paths.__repr__())
            )
    return resource_path
