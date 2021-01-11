RadioPropa
==========

RadioPropa is a fork of [CRPRopa](https://crpropa.desy.de/) to simulate radio
propagation in inhomogeneous by media ray-tracing.


# Installation
Be sure you have python, a C++ compiler and cmake!

Create a folder "RadioPropa" in the directory where you want to install it.
Create in this folder a second folder called "radiopropa". Go into this directory
using a therminal.
	$ cd [repository]/RadioPropa/radiopropa.

Clone the gitrepository in this directory
	$ git clone https://github.com/NuRadio/RadioPropa

Execute the following steps in the therminal from the radiopropa repository
	$ mkdir build
	$ cd build/
	$ cmake ..
	$ make install

Finally add the path ro the "RadioPropa" directory to your PYTHONPATH.

Done! You can now `import radiopropa` into your python scripts and use this package



# Examples
The folder `radio_example` contain some examples that can be executed from the
build directory:

 + Ice_trajectories.py `python ../radio_example/Ice.py` executes a simualtion of some ray
   trajectories starting from a given depth. Output is saved to `output_traj.h5`
 + plot_tray.py plots the data in `output_traj.h5`
 + n2linear.py simulates a single trajectory in an n2linear field.
 + plot_n2linear.py plots the field and compare the reult with the analytical
   expectation.
 + ReflectionRefraction.py `python ../radio_example/ReflectionRefraction.py`
	 executes an example with boundary layers. Output is saved to `output_traj.h5`.
