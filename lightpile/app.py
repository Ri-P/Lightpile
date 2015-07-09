# -*- coding: utf-8 -*-
"""
Usage:
    lightpile <lmc-file>
    lightpile -h | --help
    lightpile --version

Options:
    -h --help   Show this screen.
    --version   Show version.
"""
from __future__ import print_function
from __future__ import unicode_literals
import sys
import docopt

import lightpile
from lightpile.parser import Parser
from lightpile.taskcreator import Taskcreator

__author__ = "Richard Pfeifer"


def parse_commandline(argv):
    """
    Let docopt handle the CLI and return dict of arguments.
    """
    try:
        # Parse arguments
        arguments = docopt.docopt(__doc__, argv,
                                  version=lightpile.__version__)
    except docopt.DocoptExit as e:
        print(e.message)
        sys.exit(0)
    return arguments


def parse_user_lmc_file(lmc_file_path):
    """
    Parse file at unvalidated lmc_file_path to get userconf-dict.
    """
    parser = Parser()
    with open(lmc_file_path, 'r') as sf:
        userinput = sf.read()
    scanstring = userinput.decode(lightpile.codec)
    user_config_dict = parser.parse(scanstring)
    return user_config_dict


def Main():

    arguments = parse_commandline(sys.argv[1:])
    # TODO: validate CL-arguments
    lmc_file_path = arguments['<lmc-file>']

    try:
        user_configuration = parse_user_lmc_file(lmc_file_path)
        tc = Taskcreator(user_configuration)
        tasklist = tc.create_tasks()
        for task in tasklist:
            task.run()
    except KeyboardInterrupt:
        print("Shutdown requested...exiting")
    sys.exit(0)

if __name__ == "__main__":
    Main()
