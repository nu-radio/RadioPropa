import h5py
import radiopropa 

#from Ice import iceModel


#Plot Ice model
z = linspace(-300,0)

n = zeros_like(z)
dn = zeros_like(z)

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
f = h5py.File('output_traj.h5')
d = f['Trajectory3D']

figure()
s1 = subplot(111)

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
A = sqrt(d['Ax']**2 + d['Ay']**2 +d['Az']**2 )
s1.scatter(d['X'], d['Z'], c=log10(A), marker='.', s=2)
#s1.scatter(d['X'], d['Z'], c=abs(d['SN']), marker='.', s=2)
s1.plot([d['X'][0]], [d['Z'][0]], c='r', marker='*')
axhline(0., c='k')
text(1500.,5, 'Air')
text(1500.,-5, 'Ice', verticalalignment='top')
s1.set_xlabel('X [m]')
s1.set_ylabel('Z [m]')
s1.set_ylim(-300,30)
tight_layout()




#xlim(0,10)
#ylim(-240, -230)
