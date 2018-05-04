RadioPropa
==========

A fork of CRPropa to simulate radio propagation via ray tracing in inhomogeneous media.
The folder radio_example contain some examples that can be executed from the
build directory:

 + Ice_trajectories.py `python ../radio/example/Ice.py` executes a simualtion of some ray
   trajectories starting from a given depth. Output i saved to `output_traj.h5`
 + plot_tray.py plots the data in `output_traj.h5`
 + n2linear.py simulates a single trajectory in an n2linear field.
 + plot_n2linear.py plots the field and compare the reult with the analytical
   expectation.


