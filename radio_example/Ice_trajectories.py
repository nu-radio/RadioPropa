import radiopropa
import numpy as np


class ObserverZ(radiopropa.ObserverFeature):
    """
    An observer plane parallel to the xy plane in depth z.
    """
    def __init__(self, Z):
        radiopropa.ObserverFeature.__init__(self)
        self.__Z = Z

    def checkDetection(self, candidate):
        cx = candidate.current.getPosition().getZ() - self.__Z
        px = candidate.previous.getPosition().getZ() - self.__Z

        if np.sign(cx) == np.sign(px):
            return radiopropa.NOTHING
        else:
            return radiopropa.DETECTED


iceModel = radiopropa.GorhamIceModel()

if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))

    # Observer to stop imulation at =0m and z=300m
    obs = radiopropa.Observer()
    obsz = ObserverZ(0.0)
    obs.add(obsz)
    obs.setDeactivateOnDetection(True)
    sim.add(obs)

    obs2 = radiopropa.Observer()
    obsz2 = ObserverZ(-300.0)
    obs.add(obsz2)
    obs2.setDeactivateOnDetection(True)
    sim.add(obs2)

    # Output
    output = radiopropa.HDF5Output('output_traj.h5', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    output.enable(radiopropa.Output.SerialNumberColumn)
    sim.add(output)

    # Source
    source = radiopropa.Source()
    source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0, 0, -240.)))
    source.add(radiopropa.SourceParticleType(radiopropa.nucleusId(1, 1)))
    source.add(radiopropa.SourceAmplitude(1))
    source.add(radiopropa.SourceFrequency(1E6))

    #Start rays from 0 - 90 deg
    for phi in np.linspace(0,90):
        z = np.cos(phi * radiopropa.deg)
        x = np.sin(phi * radiopropa.deg)
        print('Ray Direction {} deg = ({}, 0, {})'.format(phi, x, z))
        source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
        sim.setShowProgress(True)
        sim.run(source, 1)
