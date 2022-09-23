import h5py
import radiopropa 
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#from Ice import iceModel


#Plot Ice model
z = np.linspace(-300,0)

n = np.zeros_like(z)
dn = np.zeros_like(z)

position = radiopropa.Vector3d(0,0,0) 

#for i,vz in enumerate(z):
#    position.z = vz
#    n[i] = iceModel.getValue(position)
#    dn[i] = iceModel.getGradient(position).z
#
#
#figure()
#plot(n, z)
#xlabel('Index of Refracton')
#ylabel('Z [m]')




# plot simulated data

angles = [45]
plt.figure()

for i,a in enumerate(angles):
	f = h5py.File('output_traj_'+str(a)+'_bis.h5')
	d = f['Trajectory3D']

	#SN = set(d['SN'])
	#for s in SN:
	#    idx = d['SN'] == s
	#    X = d['X'][idx]
	#    Z = d['Z'][idx]
	#    dX = X[1] - X[0]
	#    dZ = Z[1] - Z[0]
	#    phi0 = abs(arctan(dX/dZ))
	#
	#    c = cm.jet(phi0 / pi * 2)
	#    print phi0
	#
	#    s1.plot(d['X'][idx], d['Z'][idx], c=c, label='RadioPropa')
	#
	A = np.sqrt(d['Ax']**2 + d['Ay']**2 +d['Az']**2 )
	plt.scatter(d['X'], d['Z'], c=A, marker='.', s=2, norm=mpl.colors.LogNorm())
	#plt.scatter(d['X'], d['Z'], marker='.', s=2, label='Launch angle: '+str(a)+' deg')
	plt.plot([d['X'][0]], [d['Z'][0]], c='r', marker='*')
	plt.colorbar(label='Relative signal strength')
plt.axhline(0., c='k')
plt.axhline(-100., c='grey', linestyle=':')
plt.axhline(-50., c='grey', linestyle=':')
plt.axhline(-75., c='grey', linestyle=':')
plt.axhline(-15., c='grey', linestyle=':')
plt.axhline(-30., c='grey', linestyle=':')
plt.text(1500.,5, 'Air')
plt.text(1500.,-5, 'Ice', verticalalignment='top')
plt.xlabel('X [m]')
plt.ylabel('Z [m]')
plt.ylim(-300,30)
plt.tight_layout()
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(loc=4)

#plt.savefig('radio_ice_layers.png', dpi=600)
plt.show()




#xlim(0,10)
#ylim(-240, -230)
