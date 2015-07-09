# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from nose.tools import assert_equal, assert_raises

import lightpile
import lightpile.app as app
from lightpile.test.testutils import captured_out, temp_file_path

__author__ = "Richard Pfeifer"


def test_garbageargument_2():
    """
    Do SystemExit and show 'usage' for bad CLI-arguments.
    """
    argv = ['--garbage in ']
    expected_output = ("Usage:\n    lightpile <lmc-file>\n"
                       "    lightpile -h | --help\n"
                       "    lightpile --version\n")

    with captured_out() as (out, err):
        with assert_raises(SystemExit):
            app.parse_commandline(argv)
        assert_equal(expected_output, out.getvalue())


def test_versionstring_2():
    """
    CLI-call '>lightpile --version' should print correct version.
    """
    argv = ['--version']
    expected_output = lightpile.__version__
    with captured_out() as (out, err):
        with assert_raises(SystemExit):
            app.parse_commandline(argv)
        assert_equal(expected_output, out.getvalue().strip())


def test_user_dict_preparation():
    """
    Test generation of userdict from filepath of lmc-file.
    """
    content = ("[materials]\n"
               "air 1.0 0\n"
               "goo 1.7 0.1\n"
               "\xc6g\xf8 1.2 0.05\n")
    expected_dict = {
        'materials': [{'nk': (1.0+0j), 'name': "air"},
                      {'nk': (1.7+0.1j), 'name': "goo"},
                      {'nk': (1.2+0.05j), 'name': "\xc6g\xf8"}]}

    with temp_file_path(content.encode(lightpile.codec)) as filepath:
        user_conf_dict = app.parse_user_lmc_file(filepath)
        assert_equal(expected_dict, user_conf_dict)
