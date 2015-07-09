from setuptools import setup

import os.path
import re

versionfile_path = os.path.join("lightpile", "_version.py")
with open(versionfile_path, 'r') as version_file:
    version_line = version_file.read()
# ^   matches only at the beginning of the line
ver_reexpr = r"^__version__ = ['\"]([^'\"]*)['\"]"
r = re.search(ver_reexpr, version_line, re.M)
if r:
    version = r.group(1)
else:
    raise RuntimeError("Unable to find version string in {0}".format(
        versionfile_path))

config = {
    'name': 'Lightpile',
    'description': 'An optical thin film solver.',
    'author': 'Richard Pfeifer',
    'author_email': 'rip@wgd2.de',
    'url': "https://github.com/Ri-P/lightpile",
    'description': (
        "A tool to simulate dipole emitters in planar thin film stacks."),
    'license': 'GNU Affero General Public License v3 or later (AGPLv3+)',
    'long_description':
        "Lightpile is a tool to simulate dipole emitters in planar thin film"
        "stacks. This is a physics problem involving eletromagnetism and a"
        " bit of quantum mechanics. You might be interested in using Lightpile "
        "if working in the field of (organic) light emitting devices, known "
        "as OLEDs and LEDs, or quantum dots.",
    'classifiers': [
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU Affero General Public License v3 or "
        "later (AGPLv3+)",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 2.7",
        "Operating System :: POSIX :: Linux",
    ],
    'version': version,
    'install_requires': ['nose==1.3.4',
                         'numpy==1.9.2',
                         'scipy==0.15.1',
                         'matplotlib==1.4.3',
                         'pyparsing==2.0.3',
                         'docopt==0.6.1'],
    'packages': ['lightpile', 'lightpile.test'],
    'include_package_data': True,
    'entry_points': {'console_scripts': ['lightpile = lightpile.app:Main']},
    'test_suite': 'nose.collector',
}

setup(**config)
