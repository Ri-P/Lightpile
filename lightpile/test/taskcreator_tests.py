# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import print_function
import os
import os.path
from nose import with_setup

from lightpile.taskcreator import Taskcreator

result_file_paths = ["dipolestudy_data_f.txt",
                     "dipolestudy_data_p.txt",
                     "dipolestudy_graph_f.png",
                     "dipolestudy_graph_p.png"]

initial_cwd = os.getcwd()


def setup_resultfiles():
    """
    Make sure there are no previous results in directory.
    """
    are_files_in_path = [os.path.isfile(path) for path in result_file_paths]
    if True in are_files_in_path:
        raise AssertionError(
            "Abort test. There are previous results in the resultpath of this "
            "test: {0} in {1}.".format(are_files_in_path, result_file_paths))


def teardown_resultfiles():
    """
    Remove result files after test.
    """
    for path in result_file_paths:
        if os.path.isfile(path):
            os.remove(path)


def setup_change_cwd():
    lightpile_testdir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(lightpile_testdir)


def teardown_change_cwd():
    os.chdir(initial_cwd)


def setup():
    setup_resultfiles()
    setup_change_cwd()


def teardown():
    teardown_resultfiles()
    teardown_change_cwd()

fixed_angular_range = {'end': 1e4,
                       'scale': 'log',
                       'sampling_points': 50,
                       'quantity': 'u', 'begin': 1e-1,
                       'unit': 'None'}
adaptive_angular_range = {'end': 5,
                          'begin': 0,
                          'quantity': 'u',
                          'unit': 'None',
                          'sampling_points': 'adaptive',
                          'scale': None}
spectralpoint = {'unit': 'nm', 'value': 633.0, 'quantity': 'wavelength'}
ag = {'nk': (0.0767+4.366j), 'name': 'ag'}
air = {'nk': (1+0j), 'name': 'air'}
airlayer1 = {'planarmaterial': 'air', 'thickness': 0.0}
elayer_green = {'emittername': 'green'}
airlayer2 = {'planarmaterial': 'air', 'thickness': 2.5}
aglayer = {'planarmaterial': 'ag', 'thickness': 200.0}
e_green = {'iqe': 1.0, 'orientation': 'vert', 'name': 'green'}


@with_setup(setup, teardown)
def test_createtask_dipolestudy_fixed_angularrange_userconfdict():
    """
    Create a dipolestudy task with fixed angular range from a dictionary and
    let it execute.
    """
    user_conf_dict = {
        'dipolestudy': {'angularrange': fixed_angular_range.copy(),
                        'spectralpoint': spectralpoint.copy()},
        'materials': [ag.copy(), air.copy()],
        'stack': [airlayer1.copy(),
                  elayer_green.copy(),
                  airlayer2.copy(),
                  aglayer.copy()],
        'emitters': [e_green.copy()]
    }
    tc = Taskcreator(user_conf_dict)
    tasks = tc.create_tasks()
    for task in tasks:
        task.run()


@with_setup(setup, teardown)
def test_createtask_dipolestudy_adaptive_angularrange_userconfdict():
    """
    Create a dipolestudy task with adaptive angular range from a dictionary and
    let it execute.
    """
    user_conf_dict = {
        'dipolestudy': {'angularrange': adaptive_angular_range.copy(),
                        'spectralpoint': spectralpoint.copy()},
        'materials': [ag.copy(), air.copy()],
        'stack': [airlayer1.copy(),
                  elayer_green.copy(),
                  airlayer2.copy(),
                  aglayer.copy()],
        'emitters': [e_green.copy()]
    }
    tc = Taskcreator(user_conf_dict)
    tasks = tc.create_tasks()
    for task in tasks:
        task.run()
