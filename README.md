RadioPropa
==========

RadioPropa is a fork of [CRPRopa](https://crpropa.desy.de/) to simulate radio
propagation in inhomogeneous by media ray-tracing. (22/01/2018)


# Installation
Be sure you have python, a C++ compiler and cmake!

Create a folder "RadioPropa" in the directory where you want to install it.
Go into this directory using a therminal.
	$ cd [repository]/RadioPropa/

Clone the git repository in this directory
	$ git clone https://github.com/NuRadio/RadioPropa
	OR
	$ git clone git@github.com:nu-radio/RadioPropa.git

Execute the following steps in the therminal from the radiopropa repository
	$ cd radiopropa
	$ mkdir build
	$ cd build/
	$ cmake ..
	$ make install

Finally add the path to the "RadioPropa" directory to your PYTHONPATH.

Done! You can now `import radiopropa` into your python scripts and use this package



# Examples
The folder `radio_example` contain some examples that can be executed:

 + `python Ice_trajectories.py` executes a simualtion of some ray
   trajectories starting from a given depth. Output is saved to `output_traj.h5`
 + `python plot_tray.py` plots the data in `output_traj.h5`
 + `python n2linear.py` simulates a single trajectory in an n2linear field.
 + `python plot_n2linear.py` plots the field and compare the reult with the analytical
   expectation.
 + `python ReflectionRefraction.py` executes an example with boundary layers. 
   Output is saved to `output_traj.h5`.
