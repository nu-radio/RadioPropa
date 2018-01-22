import crpropa
import numpy as np
#from ObserverPlane import ObserverPlane 

kilo = 1000.



class RadioFrequency(crpropa.Module):
    """ Set the initial energy to 10 EeV """
    def __init__(self, frequency):
        crpropa.Module.__init__(self)
        self.__frequency = frequency
    def process(self, candidate):
        print 'foo'
        #candidate.setProperty("frequency", self.__frequency)
    def getDescription(self):
        return "FooBar"


# simulation setup
sim = crpropa.ModuleList()
sim.add(crpropa.SimplePropagation(1*crpropa.meter, 100*crpropa.meter))




#obs = crpropa.Observer()
#obs.add(ObserverPlane(np.asarray([0.,0., 3. * kilo*crpropa.meter]),np.asarray([10.,0, 0]), np.asarray([0,10., 0])))

obs = crpropa.Observer()
obs.add(crpropa.ObserverLargeSphere(crpropa.Vector3d(0,0,0), 99*crpropa.meter))
output = crpropa.HDF5Output('output.h5', crpropa.Output.Event3D)
output.enableProperty('frequency', 0., 'Frequency for RadioPropa')

obs.onDetection(output)
#obs.setDeactivateOnDetection(True)
sim.add(obs)

source = crpropa.Source()
source.add(crpropa.SourcePosition(crpropa.Vector3d(0, 0, 0)))
source.add(crpropa.SourceParticleType(crpropa.nucleusId(1, 1)))
source.add(crpropa.SourceEnergy(1E16 * crpropa.eV))
source.add(crpropa.SourceIsotropicEmission())


boundary = crpropa.SphericalBoundary(crpropa.Vector3d(0, 0, 0), 100*kilo*crpropa.meter)
sim.add(boundary)


sim.add(RadioFrequency(1E6))
print sim

sim.setShowProgress(True)
sim.run(source, 10000)
