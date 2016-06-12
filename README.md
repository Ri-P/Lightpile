Getting Started
===============

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

### On Linux
You can build Lightpile from source with the following major dependencies (for
a full list see requirements.txt).
- Python 2.7
- numpy (1.9.2), matplotlib (1.4.3), scipy (0.15.1)


#### Direct install with scipy available

If you have a Python interpreter with scipy installed and want to install
Lightpile into this environment it is as easy as

```bash
$ pip install git+https://github.com/ri-p/lightpile
```
I do recommend to follow the second way of using a virtual environment in
order not to clutter your python installation with numerous packages.

#### Install Lightpile into a virtual environment

It is recommended to install Lightpile into a virtual environment. Follow these
steps, assuming you start with a blank installation of Ubuntu. They have been
tested with Ubuntu 14.04 64bit.

1. Get a virtualenv with pip

    ```bash
    $ sudo apt-get install python-pip
    $ sudo pip install virtualenvwrapper
    ```
    Then add the following lines to your ~/.bashrc (for details see [here](https://virtualenvwrapper.readthedocks.org/en/laest/install.html))

    ```bash
    export WORKON_HOME=$HOME/.virtualenvs
    source /usr/local/bin/virtualenvwrapper.sh
    ```
    and create your virtual environment

    ```bash
    $ source ~/.bashrc
    $ mkvirtualenv env-lp
    ```

2. Install numpy, scipy, matplotlib

    ```bash
    (env-lp)$ sudo apt-get install python-dev g++ libpng-dev libfreetype6-dev libjpeg8-dev liblapack-dev gfortran
    (env-lp)$ pip install numpy==1.9.2
    (env-lp)$ pip install matplotlib==1.4.3
    (env-lp)$ pip install scipy==0.15.1
    ```

3. Install Lightpile

    ```bash
    (env-lp)$ sudo apt-get install git
    (env-lp)$ pip install git+https://github.com/ri-p/lightpile
    ```

### On Windows
Lightpile depends on the numpy-package. Due to the difficulties of compiling numpy from source on Windows we suggest using a Python distribution with precompiled scientific packages, i.e. [Anaconda] (https://www.continuum.io/downloads) or [Canopy] (https://www.enthought.com/products/canopy). The following instruction uses Anaconda.

1. Install Anaconda and prepare a new conda-environment 'env-lp'

	Download [Python 2.7 Anaconda for Windows] (https://www.continuum.io/downloads) and install.

	```bash
	$ conda create --name env-lp python==2.7 numpy==1.9.2 scipy==0.15.1 matplotlib==1.4.3
	```

2. Install Lightpile

	```bash
	$ activate env-lp
	(env-lp)$ pip install https://github.com/ri-p/lightpile/archive/master.zip
	```


Minimal Example
---------------

As an introduction we reproduce the situation shown in Fig. 10 of the paper of
[Ford and Weber](http://deepblue.lib.umich.edu/bitstream/handle/2027.42/24649/0000062.pdf)

In your favourite text editor create a file with the following content and name
it 'ford_fig10.lps' (lightpile scene)
```
[materials]
# name          n       k
air             1.      0.
dielectric      1.4142  0
ag              0.0767  4.366

[emitters]
# name          orientation     iqe
red_emitter     vert            1.0

[stack]
air             0
red_emitter
air             0.2
dielectric      800
ag              0

[dipolestudy]
spectralpoint wavelength nm 633
angularrange u None 0.1 1.6 adaptive

[graph_p_a]
yscale log
ylim 1e-3 1e3
```

Let's start the calculation
```
(env-lp)$ lightpile ford_fig10.lps
```

Lightpile creates for output files
* dipolestudy_graph_f.png
* dipolestudy_graph_p.png
* dipolestudy_data_f.txt
* dipolestudy_data_p.txt

Compare dipolestudy_graph_p.png with figure 10 of the paper.

For explanations concerning the physical meaning of p(u) and f(u) please read
the detailed [Lightpile examples](https://ri-p.github.io/Lightpile/examples.html).

Once you are done with your calculations, leave the environment with
```
(env-lp)$ deactivate
```
