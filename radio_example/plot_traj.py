import h5py
from mpl_toolkits.mplot3d import axes3d
from radiopropa import Mpc

f = h5py.File('output_traj.h5')
d = f['Trajectory3D']


N = min(d.size, 10000)

fig = plt.figure()
ax = fig.gca(projection='3d')# , aspect='equal'
ax.scatter(d[:N]['X'] * Mpc,d[:N]['Y'] * Mpc, d[:N]['Z'] * Mpc, 'o', c=d[:N]['E'] )

R = 10
ax.set_zlim(-R, R)
ax.set_ylim(-R, R)
ax.set_zlim(-6, 6)
