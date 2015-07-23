Getting Started
===============

Overview
--------

Lightpile is a tool to simulate dipole emitters in thin film stacks.
This is a physics problem involving eletromagnetism and a bit of quantum
mechanics. You might be interested in using Lightpile if working in the field
of (organic) light emitting devices, known as (O)LEDs, or quantum dots.

At the moment Lightpile calculates the emission behaviour of a monochromatic
dipole emitter of arbitrary orientation within a pile of thin films.

It is written by Richard Pfeifer and licensed under the GNU AGPL license.

Installation
------------

### On Linux
You can build Lightpile from source with the following major dependencies
- Python 2.7
- numpy, matplotlib, scipy

For a full list see requirements.txt.

Build and install Lightpile with

`pip install git+git://github.com/ri-p/lightpile`

#### Install Lightpile into a virtual environment
It is recommended to install Lightpile into a virtual environment.
The following steps have been tested with Ubuntu 14.04 64bit

1. Get a virtualenv with pip
  ```
  $ sudo apt-get install python-pip
  $ sudo pip install virtualenvwrapper
  ```
  Add the following lines to your ~/.bashrc (for details see
[here](https://virtualenvwrapper.readthedocks.org/en/laest/install.html) ):
  ```
  export WORKON_HOME=$HOME/.virtualenvs
  source /usr/local/bin/virtualenvwrapper.sh
  ```
  and create your virtual environment
  ```
  $ source ~/.bashrc
  $ mkvirtualenv env-lp
  ```
2. Install numpy, scipy, matplotlib
  ```
  (env-lp)$ sudo apt-get install python-dev g++ libpng-dev libfreetype6-dev libjpeg8-dev liblapack-dev gfortran
  (env-lp)$ pip install numpy
  (env-lp)$ pip install matplotlib
  (env-lp)$ pip install scipy
  ```
3. Install Lightpile
  ```
  (env-lp)$ sudo apt-get install git
  (env-lp)$ pip install git+https://github.com/ri-p/lightpile
  ```

### On Windows
We plan to release binaries for Windows.


Minimal Example
---------------

As an introduction we reproduce the situation shown in Fig. 10 of the paper of
[Ford and Weber](http://deepblue.lib.umich.edu/bitstream/handle/2027.42/24649/0000062.pdf)

In your favourite text editor create a file with the follwing content and name
it 'ford_fig10.lps' (lightpile scene)
```
[materials]
# name		n	k
air		1. 	0.
dielectric	1.4142	0
ag		0.0767	4.366

[emitters]
# name		orientation 	iqe
red_emitter	vert		1.0

[stack]
air		0
red_emitter
air		0.2
dielectric	800
ag		0

[dipolestudy]
spectralpoint wavelength nm 633
angularrange u None 0.1 1.6 adaptive
```

Let's start the calculation
```
(env-lp)$ lightpile ford_fig10.lps

Lightpile creates for output files
* dipolestudy_graph_f.png
* dipolestudy_graph_p.png
* dipolestudy_data_f.txt
* dipolestudy_data_p.txt

For explanation concerning the physical meaning of p(u) and f(u) please read
the detailed [Lightpile examples](https://ri-p.github.io/Lightpile/examles.html).
