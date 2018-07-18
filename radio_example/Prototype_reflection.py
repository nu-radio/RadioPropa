import radiopropa
import numpy as np


def printVector(v):
    print('v = {}, {}, {}; |v| = {}'.format(v.x, v.y, v.z, v.getR()))
    if v.x != v.x or v.y != v.y or v.z != v.z:
        brk
    if v.getR() > 2:
        brk


class Discontinuity(radiopropa.Module):
    """
    """
    def __init__(self, surface, n1, n2):
        radiopropa.Module.__init__(self)
        self.__surface = surface
        self.__n1 = n1
        self.__n2 = n2

    def process(self, candidate):

        cx = self.__surface.distance(candidate.current.getPosition())
        px = self.__surface.distance(candidate.previous.getPosition())

        if np.sign(cx) == np.sign(px):
            candidate.limitNextStep(abs(cx))
        else:
            ## The secondary propagates further, while the candidate is
            ## reflected: legacy from CRPropa interface as secondaries have same
            ## direction as parents.

            # calculate inersection point
            dp = candidate.current.getPosition() -  candidate.previous.getPosition()
            intersectionPoint = candidate.previous.getPosition() + dp.getUnitVector() * px

            #surface normal in intersection point 
            localNormal= self.__surface.normal(intersectionPoint)

            # angle to the surface normal
            #alpha = candidate.current.getDirection().getAngleTo(localNormal)
            ca = abs(candidate.current.getDirection().dot(localNormal))
            alpha = np.arccos(ca)

            if px < 0:
                 N1 = self.__n1
                 N2 = self.__n2
            else:
                 N1 = self.__n2
                 N2 = self.__n1

            # reflection according to Snell's law
            #ca = np.cos(alpha)
            sa = abs(np.sin(alpha))
            sb = N1 / N2 * sa

            # vector in plane ip perpendicular to direction and normal
            ip = localNormal.cross(candidate.current.getDirection())
            A = candidate.current.getAmplitude()
            Aperp = localNormal * A.dot(localNormal) 
            Apara = ip * A.dot(ip)
            Aperp = A - Apara

            print('  ')
            print('  Aperp:')
            printVector(Aperp)
            print('  Apara:')
            printVector(Apara)
            print('  Sum')
            printVector(Apara + Aperp)
            print('  ')


            if abs(sb) <1:
                # No total reflection, calculate reflection and transmission
                # coefficient from Fresnell equations

                cb = abs(np.cos(np.arcsin(sb)))

                T_perp = 2 * N1 * ca / (N1 * ca + N2 * cb)
                R_perp = (N1*ca - N2*cb) / (N1 * ca + N2 * cb)

                T_para = 2 * N1 * ca / (N2 * ca + N1 * cb)
                R_para = (N2 * ca - N1 * cb) / (N2 * ca + N1 * cb)
                
                print('Fresnell:')
                print('   ca = {}, cb = {}'.format(ca, cb))
                print('  T_perp = {}, T_para = {}'.format(T_perp, T_para))
                print('  R_perp = {}, R_para = {}'.format(R_perp, R_para))

                c2 = candidate.clone()
                Vnew = c2.current.getDirection()
                phi = Vnew.getPhi()
                Vnew.setRThetaPhi(1., np.arcsin(sb), phi)
                c2.current.setDirection(Vnew)
                beta = np.arcsin(sb)
                Aperp_p = Aperp.getRotated(ip, beta-alpha)
                TransmittedAmplitude = Apara * T_para + Aperp_p * T_perp  
                # correction factor to account for the increased beamwidth
                c  = np.sqrt(abs(N2 / N1 * cb/ca))
                print("N2 = {} N1 = {} cb = {}  ca = {} c = {} ".format(N2, N1, cb, ca, c))
                print('Transmission:')
                printVector(TransmittedAmplitude)
                TransmittedAmplitude *= c
                printVector(TransmittedAmplitude)

                c2.current.setAmplitude(TransmittedAmplitude)
                candidate.addSecondary(c2)
            else:
                # Total reflection, no transmission
                R_perp = 1.
                R_para = 1.



            # resnell equations
            #sb = sel.__n1 / self.__n2 * sa
            #b = arcsin(sb)

            ## perpendicular r,t coefficients
            #rper = -1. * sin(a-b) / sin(a+b) 
            #tper = 2* sb*ca / sin(a+b)
            ## parallel r,t coefficients
            #rpar = tan(a-b) / tan(a+b) 
            #tpar = 2*sb*sa / sin(a+b) / cos(a-b)


            ## reflection
            V = candidate.current.getDirection()
            u = localNormal * (V.dot(localNormal))
            new_direction = V - u*2
            candidate.current.setDirection(new_direction)

            print ('Rpara: {}, Rperp: {}'.format(R_para, R_perp))

            Aperp_p = Aperp - localNormal * (Aperp.dot(localNormal)) * 2
            ReflectedAmplitude = Apara * R_para + Aperp_p * R_perp

            #ReflectedAmplitude= radiopropa.Vector3d(1,0,0)
            candidate.current.setAmplitude(ReflectedAmplitude)
            print('Reflected Amplitude:')
            printVector(ReflectedAmplitude)

            # update position slightly to move on correct side of plane
            X = candidate.current.getPosition()
            candidate.current.setPosition(X + new_direction * candidate.getCurrentStep())



iceModel = radiopropa.GorhamIceModel(z0=-10.)
airBoundary = Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,0), radiopropa.Vector3d(0,0,1)), 1.5 , 1.)
firnLayer = Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-10.), radiopropa.Vector3d(0,0,1)), iceModel.getValue(radiopropa.Vector3d(0,0,-10.001)) , 1.5)


alpha = 30.

c = radiopropa.Candidate()
c.previous.setPosition(radiopropa.Vector3d(0,0,-.1))
c.previous.setDirection(radiopropa.Vector3d(0,0,1))
c.current.setDirection(radiopropa.Vector3d(0,0,1))
amplitude = radiopropa.Vector3d(1,1,0).getUnitVector()
c.previous.setAmplitude(amplitude)
c.current.setAmplitude(amplitude)
print('Start')
printVector(c.current.getDirection())
printVector(c.current.getAmplitude())
#print(c.current.getAmplitude().getAngleTo(c.current.getDirection()))

direction = radiopropa.Vector3d()
direction.setRThetaPhi(1, alpha/180.*np.pi, 0)
c.previous.setDirection(direction)
c.current.setDirection(direction)
c.current.setPosition(c.previous.getPosition() + direction)

print('')
print('Before:')
printVector(c.current.getDirection())
printVector(c.current.getAmplitude())
print(c.current.getAmplitude().getAngleTo(c.current.getDirection()))

airBoundary.process(c)
print('')
print('After:')
print("refl: ")
printVector(c.current.getDirection())
printVector(c.current.getAmplitude())
print(c.current.getAmplitude().getAngleTo(c.current.getDirection()))

if len(c.secondaries) > 0:
    print("trans: ")
    printVector(c.secondaries[0].current.getDirection())
    printVector(c.secondaries[0].current.getAmplitude())
    print(c.secondaries[0].current.getAmplitude().getAngleTo(c.secondaries[0].current.getDirection()))




if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))

    sim.add(airBoundary)
    sim.add(firnLayer)

    sim.add(radiopropa.MinimumAmplitude(.05))
    sim.add(radiopropa.MaximumTrajectoryLength(1500 * radiopropa.meter))

    # Observer to stop imulation at = 30m and z=300m
    obs = radiopropa.Observer()
    obsz = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,30), radiopropa.Vector3d(0,0,1)))
    obs.add(obsz)
    obs.setDeactivateOnDetection(True)
    sim.add(obs)

    obs2 = radiopropa.Observer()
    obsz2 = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,-300), radiopropa.Vector3d(0,0,1)))
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
    source.add(radiopropa.SourceAmplitude(1))
    source.add(radiopropa.SourceFrequency(1E6))

    #Start rays from 0 - 90 deg
    for phi in np.linspace(0,90, 50):
        z = np.cos(phi * radiopropa.deg)
        x = np.sin(phi * radiopropa.deg)
        print('Ray Direction {} deg = ({}, 0, {})'.format(phi, x, z))
        source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
        sim.setShowProgress(True)
        sim.run(source, 1)
