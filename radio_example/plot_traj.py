import h5py
from mpl_toolkits.mplot3d import axes3d
from radiopropa import Mpc

f = h5py.File('output_traj.h5')
d = f['Trajectory3D']


N = min(d.size, 10000)
dN = 1 #d.size / N
idx = arange(0,d.size, dN)

#fig = plt.figure()
#ax = fig.gca(projection='3d')# , aspect='equal'
#ax.scatter(d[::dN]['X'],d[::dN]['Y'], d[::dN]['Z'], 'o', c=d[::dN]['D'] )
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')


figure()
plot(d['X'], d['Z'], 'ro')
#xlim(0,10)
#ylim(-240, -230)
