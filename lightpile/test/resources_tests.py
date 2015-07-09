# -*- coding: utf-8 -*-

import os
import os.path
from nose import with_setup
from nose.tools import assert_raises, assert_equal

from lightpile.resources import get_resource_path

initial_cwd = os.getcwd()


def setup_cwd():
    # set cwd to lightpile test dir
    lightpile_testdir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(lightpile_testdir)


def teardown_cwd():
    # set cwd back
    os.chdir(initial_cwd)


@with_setup(setup_cwd, teardown_cwd)
def test_relative_path():

    magicstring = '888777666'
    rsc_relative_path = "data/green.pl"
    resource_path = get_resource_path(rsc_relative_path)
    assert(os.path.isabs(resource_path))

    # test for the right file using magic number
    with open(resource_path, 'r') as rsc_file:
        magicstring_read = rsc_file.readline().split()[4]
    assert_equal(magicstring_read, magicstring)


def test_invalid_path():
    rsc_invalid_path = 543
    assert_raises(Exception, get_resource_path, rsc_invalid_path)
