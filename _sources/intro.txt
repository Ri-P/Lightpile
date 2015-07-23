
Getting Started
===============

.. contents::
   :local:

Overview
--------

Lightpile is a tool to simulate dipole emitters in thin film stacks.
This is a physics problem involving eletromagnetism and a bit of quantum
mechanics. You might be interested in using Lightpile if working in the field
of (organic) light emitting devices, known as (O)LEDs, or quantum dots.

At the moment Lightpile calculates the emission behaviour of a monochromatic
dipole emitter of arbitrary orientation within a stack (pile) of thin films.

It is written by Richard Pfeifer and licensed under the GNU AGPL license.


Installation
------------

On Linux
^^^^^^^^
You can build Lightpile from source with the following major dependencies (for
a full list see requirements.txt)

* Python 2.7
* numpy (1.9.2), matplotlib (1.4.3), scipy (0.15.1)

Direct install with scipy available
"""""""""""""""""""""""""""""""""""

If you have a Python interpreter with scipy installed and want to install
Lightpile into this environment it is as easy as

.. code-block:: bash

  $ pip install git+https://github.com/ri-p/lightpile

I do recommend to follow the second way of using a virtual environment in
order not to clutter your python installation with numerous packages.

Install Lightpile into a virtual environment
""""""""""""""""""""""""""""""""""""""""""""

It is recommended to install Lightpile into a virtual environment. Follow these
steps, assuming you start with a blank installation of Ubuntu. They have been
tested with Ubuntu 14.04 64bit.

1. Get a virtualenv with pip

  .. code-block:: bash

    $ sudo apt-get install python-pip
    $ sudo pip install virtualenvwrapper

  Then add the following lines to your ~/.bashrc (for details see `here`_)

  .. _here: https://virtualenvwrapper.readthedocks.org/en/laest/install.html

  .. code-block:: bash

     export WORKON_HOME=$HOME/.virtualenvs
     source /usr/local/bin/virtualenvwrapper.sh

  and create your virtual environment

  .. code-block:: bash

    $ source ~/.bashrc
    $ mkvirtualenv env-lp

2. Install numpy, scipy, matplotlib

  .. code-block:: bash

     (env-lp)$ sudo apt-get install python-dev g++ libpng-dev libfreetype6-dev libjpeg8-dev liblapack-dev gfortran
     (env-lp)$ pip install numpy==1.9.2
     (env-lp)$ pip install matplotlib==1.4.3
     (env-lp)$ pip install scipy==0.15.1

3. Install Lightpile

  .. code-block:: bash

     (env-lp)$ sudo apt-get install git
     (env-lp)$ pip install git+https://github.com/ri-p/lightpile

On Windows
^^^^^^^^^^

We plan to release binaries for Windows.


Minimal Example
---------------

As an introduction we reproduce the situation shown in Fig. 10 of the paper of
`Ford and Weber`_

.. _Ford and Weber: http://deepblue.lib.umich.edu/bitstream/handle/2027.42/24649/0000062.pdf

In your favourite text editor create a file with the follwing content and name
it 'ford_fig10.lps' (lightpile scene)::

    # Vertical dipole above dielectric layer on silver mirror.
    #
    [materials]
    # name	n	k
    air		1. 	0.
    dielectric	1.4142	0.
    ag		0.0767	4.366

    [emitters]
    # name	orientation 	iqe
    red_emitter	vert		1.0

    [stack]
    # describe stack using 'materials' and 'emitters' names
    air		0       # materialname  thickness in nm
    red_emitter         # emittername
    air		0.2
    dielectric	800
    ag		0

    [dipolestudy]
    # we start a task 'dipolestudy'. It needs to arguments:
    # 1) the wavelength of the emitter is 633nm 
    spectralpoint wavelength nm 633
    # 2) the intervall of directions or modes is: u = 0.1 ... 1.6 (no unit) 
    #    with an automatic adaptive sampling point spacing
    angularrange u None 0.1 1.6 adaptive

Let's start the calculation

.. code-block:: bash

   (env-lp)$ lightpile ford_fig10.lps


Lightpile creates for output files

* dipolestudy_graph_f.png
* dipolestudy_graph_p.png
* dipolestudy_data_f.txt
* dipolestudy_data_p.txt

For explanation concerning the physical meaning of p(u) and f(u) please read
the detailed `Lightpile examples`_.

.. _Lightpile examples: https://ri-p.github.io/Lightpile/examples.html
