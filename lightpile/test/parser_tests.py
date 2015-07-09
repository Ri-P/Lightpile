# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from nose.tools import assert_equal, assert_raises

import lightpile.parser as parser
from lightpile.test.testutils import captured_err

# Test:
# - only sections but with UNIX and Windows lineending
#


class ParseSample(object):

    def __init__(self, lmc_string, expected_dict):
        self.string = lmc_string
        self.exp_dict = expected_dict

lmc_str_dipolestudy_ford = """\
[materials]
# name  n       k
ag      0.0767  4.366
air     1.

[emitters]   # define emitters here
# name  orientation iqe
green   vert        1.0

[stack]
# define stack here
air     0
green
air     0.05
ag      200

[dipolestudy]
spectralpoint wavelength nm 633
angularrange u None 0.1 1e4 1000 log

[graph_p_a]
yscale log
lines all
xlim 0.1 1e4
ylim 1e-2 1e4
externaldata "examples/ex_ford1984/ford1984_fig4a.dat" "Ford d=0.05nm"

[graph_f_a]
ylim auto auto
xlim 0.1 1e4
xscale log
lines all
externaldata "examples/ex_ford1984/ford1984_fig4a.dat" "Ford d=0.05nm"
"""

expected_dict_dipolestudy_ford = {
    'dipolestudy': {'angularrange': {'end': 1e4, 'scale': 'log',
                                     'sampling_points': 1000,
                                     'quantity': 'u', 'unit': 'None',
                                     'begin': 0.1},
                    'spectralpoint': {'unit': 'nm',
                                      'quantity': 'wavelength',
                                      'value': 633}},
    'materials': [{'nk': (0.0767+4.366j), 'name': 'ag'},
                  {'nk': (1+0j), 'name': 'air'}],
    'emitters': [{'iqe': 1.0,
                  'orientation': 'vert',
                  'name': 'green'}],
    'stack': [{'planarmaterial': 'air', 'thickness': 0.0},
              {'emittername': 'green'},
              {'planarmaterial': 'air', 'thickness': 0.05},
              {'planarmaterial': 'ag', 'thickness': 200.}],
    'graph_p_a': {'yscale': 'log',
                  'lines': 'all',
                  'xlim_left': 0.1,
                  'xlim_right': 10000.,
                  'ylim_bottom': 0.01,
                  'ylim_top': 10000,
                  'externaldata': {
                      'path': 'examples/ex_ford1984/ford1984_fig4a.dat',
                      'label': 'Ford d=0.05nm'}},
    'graph_f_a': {'xscale': 'log',
                  'lines': 'all',
                  'xlim_left': 0.1,
                  'xlim_right': 10000,
                  'ylim_bottom': 'auto',
                  'ylim_top': 'auto',
                  'externaldata': {
                      'path': 'examples/ex_ford1984/ford1984_fig4a.dat',
                      'label': 'Ford d=0.05nm'}},
}

parsesample_dipolestudy_ford = ParseSample(
    lmc_str_dipolestudy_ford,
    expected_dict_dipolestudy_ford)

parsesample_material_specialcharacters = ParseSample(
    """\
    [materials]
    grünerSchleiµ 1.4 5.6
    große_moleküle 2.2
    CH5:NDP+ 1.7 0.001
    µüsomorph-3%5 0.9
    """,
    {'materials': [{'nk': (1.4+5.6j), 'name': 'grünerSchleiµ'},
                   {'nk': (2.2+0j), 'name': 'große_moleküle'},
                   {'nk': (1.7+0.001j), 'name': 'CH5:NDP+'},
                   {'nk': (0.9+0j), 'name': 'µüsomorph-3%5'}]})


def test_parser_fail_on_nonunicode_string():
    """
    Parser should only accept unicodestring and fail on bytestring.
    """
    p = parser.Parser()
    bytestring = str("[materials]\nmeinMaterial 1.4 5.6\nMaterial2 1.7 0.01\n")
    with assert_raises(TypeError):
        p.parse(bytestring)


def test_fail_on_bad_input():
    """
    Lightpile should exit on bad input with sane error message on stderr.
    """
    p = parser.Parser()
    badinput = "[materials]\nmat1 2 3 4"
    expected_stderr = ('1\tmat1 2 3 4\nExpected end of line (at char 9), '
                       '(line:1, col:10)\nPossible cause: Unknown or mistyped '
                       'keyword after given expected end of text.\n')
    with captured_err() as err:
        with assert_raises(SystemExit) as cm:
            p.parse(badinput)
            assert_equal(cm.exception.code, 1)
        assert_equal(err.getvalue(), expected_stderr)


def test_parser_dipolestudy_ford():
    """
    Test parsing a complete example lmc-file.
    """
    p = parser.Parser()
    tokendict = p.parse(parsesample_dipolestudy_ford.string)
    for key in tokendict.iterkeys():
        assert_equal(tokendict[key],
                     parsesample_dipolestudy_ford.exp_dict[key])


def test_parser_specialcharacters():
    """
    Test some non-ascii unicode characters as identifiers in parse string.
    """
    p = parser.Parser()
    tokendict = p.parse(parsesample_material_specialcharacters.string)
    for key in tokendict.iterkeys():
        assert_equal(tokendict[key],
                     parsesample_material_specialcharacters.exp_dict[key])
