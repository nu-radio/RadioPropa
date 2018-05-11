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


class DiscontinuityZ(radiopropa.Module):
    """
    A layer defined by point X0 and vectors V1, V2. Part of a candidate
    crossing the layer is reflected. 
    Future propagation in scalar field will take care of directivity change automatically???
    """
    def __init__(self, Z, n1, n2):
        radiopropa.Module.__init__(self)
        self.__Z = Z
        self.__n1 = n1
        self.__n2 = n2

    def process(self, candidate):

        cx = candidate.current.getPosition().getZ() - self.__Z
        px = candidate.previous.getPosition().getZ() - self.__Z

        if np.sign(cx) == np.sign(px):
            candidate.limitNextStep(abs(cx))
        else:
            #E = candidate.current.getAmplitude()

            ## The secondary propagates further, while the candidate is
            ## reflected: legacy from CRPropa interface as secondaries have same
            ## direction as parents.

            #snells law
            if px < 0:
                sin_thetaNew = self.__n1 / self.__n2 * np.sin(candidate.current.getDirection().getTheta())
            else:
                sin_thetaNew = self.__n2 / self.__n1 * np.sin(candidate.current.getDirection().getTheta())
            
            if abs(sin_thetaNew) <= 1:
                # No total reflection
                c2 = candidate.clone()
                Vnew = c2.current.getDirection()
                phi = Vnew.getPhi()
                Vnew.setRThetaPhi(1., np.arcsin(sin_thetaNew), phi)
                c2.current.setDirection(Vnew)
                candidate.addSecondary(c2)

            V = candidate.current.getDirection()

            # resnell equations
            #ca = V.getZ()
            #a = V.getTheta()
            #sa = sin(a)
            #sb = sel.__n1 / self.__n2 * sa
            #b = arcsin(sb)

            ## perpendicular r,t coefficients
            #rper = -1. * sin(a-b) / sin(a+b) 
            #tper = 2* sb*ca / sin(a+b)
            ## parallel r,t coefficients
            #rpar = tan(a-b) / tan(a+b) 
            #tpar = 2*sb*sa / sin(a+b) / cos(a-b)


            # reflection
            self.__normal = radiopropa.Vector3d(0,0,1)
            u = self.__normal * (V.dot(self.__normal))
            new_direction = V - u*2
            candidate.current.setDirection(new_direction)

            # update position slightly to move on correct side of plane
            X = candidate.current.getPosition()

            candidate.current.setPosition(X + new_direction * candidate.getCurrentStep())
            # update position (this is a hack to avoid double scatter)
            #candidate.previous.setPosition(candidate.current.getPosition())






iceModel = radiopropa.GorhamIceModel()

if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))

    iceAirBoundary = DiscontinuityZ(0., iceModel.getValue(radiopropa.Vector3d(0,0,-1E-128)), 1.)
    sim.add(iceAirBoundary)

    # Observer to stop imulation at =0m and z=300m
    obs = radiopropa.Observer()
    obsz = ObserverZ(30.0)
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
