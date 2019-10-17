import radiopropa
#from ObserverPlane import ObserverPlane 

class MyTestMaterial(radiopropa.ScalarField):
    """
    No physics in here, just a demo.
    """
    def __init__(self):
        radiopropa.ScalarField.__init__(self)
    def add(self, volume, field):
        self.components.append([volume, field])

    def getValue(self, position):
        print 'foo'
        return 2 + np.sin(position.x)

    def getGradient(self, position):
        print 'bar'
        return radiopropa.Vector3d(np.cos(position.x), 0, 0)




# simulation setup
sim = radiopropa.ModuleList()

field = MyTestMaterial() # It is important to declear python objects once in the code, as there will be a segfault if done directly in the constructor of another object. Stupid swig bug. 
sim.add(radiopropa.PropagationCK(field, 1E-20, .0001, 1.))




#obs = radiopropa.Observer()
#obs.add(ObserverPlane(np.asarray([0.,0., 3. * kilo*radiopropa.meter]),np.asarray([10.,0, 0]), np.asarray([0,10., 0])))

#obs = radiopropa.Observer()
#obs.add(radiopropa.ObserverLargeSphere(radiopropa.Vector3d(0,0,0), 99*radiopropa.meter))
output = radiopropa.HDF5Output('output_traj.h5', radiopropa.Output.Trajectory3D)
output.setLengthScale(radiopropa.meter)
#output.enableProperty('frequency', 0., 'Frequency for RadioPropa')

#obs.onDetection(output)
#obs.setDeactivateOnDetection(True)
sim.add(output)

source = radiopropa.Source()

source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0, 0, 0)))
source.add(radiopropa.SourceAmplitude(1))


#source.add(radiopropa.SourceIsotropicEmission())
source.add(radiopropa.SourceDirection(radiopropa.Vector3d(.3,0,1.)))


# Not constructiong this outside the add method will cause segfault
source.add(radiopropa.SourceFrequency(1E6))

#Two transmissive layers at +/- 10 m

sim.add(radiopropa.MinimumAmplitude(1E-2))
sim.add(radiopropa.MaximumTrajectoryLength(2000))

#boundary = radiopropa.SphericalBoundary(radiopropa.Vector3d(0, 0, 0), 100*radiopropa.meter)
#sim.add(boundary)



sim.setShowProgress(True)
sim.run(source, 1)
#print rf
