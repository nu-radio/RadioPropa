import radiopropa
import numpy as np
#from ObserverPlane import ObserverPlane 

kilo = 1000.



class RadioFrequency(radiopropa.Module):
    """ Set the initial energy to 10 EeV """
    def __init__(self, frequency):
        radiopropa.Module.__init__(self)
        self.__frequency = frequency
    def process(self, candidate):
        print 'foo'
        #candidate.setProperty("frequency", self.__frequency)
    def getDescription(self):
        return "FooBar"


# simulation setup
sim = radiopropa.ModuleList()
sim.add(radiopropa.SimplePropagation(1*radiopropa.meter, 100*radiopropa.meter))




#obs = radiopropa.Observer()
#obs.add(ObserverPlane(np.asarray([0.,0., 3. * kilo*radiopropa.meter]),np.asarray([10.,0, 0]), np.asarray([0,10., 0])))

obs = radiopropa.Observer()
obs.add(radiopropa.ObserverLargeSphere(radiopropa.Vector3d(0,0,0), 99*radiopropa.meter))
output = radiopropa.HDF5Output('output.h5', radiopropa.Output.Event3D)
output.enableProperty('frequency', 0., 'Frequency for RadioPropa')

obs.onDetection(output)
#obs.setDeactivateOnDetection(True)
sim.add(obs)

source = radiopropa.Source()
source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0, 0, 0)))
source.add(radiopropa.SourceParticleType(radiopropa.nucleusId(1, 1)))
source.add(radiopropa.SourceEnergy(1E16 * radiopropa.eV))
source.add(radiopropa.SourceIsotropicEmission())


boundary = radiopropa.SphericalBoundary(radiopropa.Vector3d(0, 0, 0), 100*kilo*radiopropa.meter)
sim.add(boundary)

# Not constructiong this outside the add method will cause segfault
#rf = RadioFrequency(1E6)
#sim.add(rf)
print sim

sim.setShowProgress(True)
sim.run(source, 10000)
#print rf
