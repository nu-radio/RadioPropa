import h5py
import radiopropa 

from n2linear import iceModel, n0, a


figure()
s1 = subplot(221)

# plot simulated data
f = h5py.File('output_traj.h5')
d = f['Trajectory3D']
s1.plot(d['X'], d['Z'], 'ro', label='RadioPropa')

# plot analytical expectation
phi0 = 1. * arctan(1.)
xi0s = n0 * cos(phi0)

z = d['Z'] 
x = 2 *xi0s / a * (sqrt(-xi0s**2+n0**2 +a*z) - sqrt(-xi0s**2 + n0**2) )

s1.plot(x, z, label='Analytical')
s1.legend(loc='upper left', fontsize='small')


# plot difference to analytical expectation
s3 = subplot(223)
s3.plot(x, abs((x-d['X'])) / x)

s3.set_xlabel('X [m]')

s1.set_ylabel('Z [m]')
s3.set_ylabel('$\\Delta$ X / X')
s1.set_xticklabels([])
s3.semilogy()

tight_layout()
subplots_adjust(hspace=0.08, wspace=0.08)


# plot refractivity profile`
n = zeros_like(z)
dn = zeros_like(z)
position = radiopropa.Vector3d(0,0,0) 

for i,vz in enumerate(z):
    position.z = vz
    n[i] = iceModel.getValue(position)
    dn[i] = iceModel.getGradient(position).z

s2 = subplot(222)
s2.plot(n, z)
s2.set_yticklabels([])
s2.set_xlabel('Refractive Index')
