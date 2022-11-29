RadioPropa
==========

RadioPropa is a fork of [CRPRopa](https://crpropa.desy.de/) to simulate radio
propagation in inhomogeneous by media ray-tracing. (22/01/2018)


# Installation
## General
Be sure you have python, a C++ compiler with c++11 support (gcc, clang and icc are known to work), swig (3.0.4 or higher), cmake and HDF5!

Go into the directory where you want to install it using a therminal.

        $ cd [software_repository]

Clone the git repository in this directory

        $ git clone https://github.com/NuRadio/RadioPropa.git
        OR
        $ git clone git@github.com:nu-radio/RadioPropa.git

Execute the following steps in the therminal from the radiopropa repository

        $ cd RadioPropa
        $ cd radiopropa
        $ mkdir build
        $ cd build/
        $ cmake ..
        $ make install

Finally add the path to the "RadioPropa" directory `[software_repository]/RadioPropa` to your PYTHONPATH.

Done! You can now `import radiopropa` into your python scripts and use this package

## Resolving issues
The installation may sometimes require some extra arguments or steps depending on your local installation. Here you will find some common issues:

+ For Apple users, it is possible that the compiler cannot find you python installation. You can resolve this by handing the location to your local python installation to the cmake command by adding some flags as shown here:

       $ cmake -DPYTHON_LIBRARY=$(python3-config --prefix)[path to 'libpython3.10.dylib'] 
        -DPYTHON_INCLUDE_DIR=$(python3-config --prefix)[path to '/include/python3.10'] ..

+ If you installed python via homebrew, you can adjust cmake as follows:

       $ cmake -DPYTHON_LIBRARY=/opt/homebrew/Cellar/python@3.10/3.10.8/Frameworks/Python.framework/Versions/3.10/lib/libpython3.10.dylib  -DPYTHON_INCLUDE_DIR=/opt/homebrew/Cellar/python@3.10/3.10.8/bin/python3.10 ..

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
