import radiopropa
import numpy as np
import radiotools.helper as hp

class DictOutput(radiopropa.Output):
    def __init__(self,candidate,dict,trajectory3D=True):
        _dict = dict
        _candidate = candidate
        if trajectory3D:
            _dict['path_3D'] = np.array([])

    def process(self,candidate):
        _dict['path_3D'].append(candidate.getProperty('pos'))



airBoundary = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,0), radiopropa.Vector3d(0,0,1)), 1.3, 1)

iceModel = radiopropa.GorhamIceModel()

if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))

    # Observer to stop imulation at =0m and z=300m
    sim.add(airBoundary)

    obs2 = radiopropa.Observer()
    obsz2 = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,-300), radiopropa.Vector3d(0,0,1)))
    obs2.add(obsz2)
    obs2.setDeactivateOnDetection(True)
    sim.add(obs2)



    # Source
    x1 = [0, 0, -240]
    x1_dir = (80,0)
    channel = [100, 0, -100]
    
    phi_direct, theta = hp.cartesian_to_spherical(*(np.array(channel)-np.array(x1)))
    phi_direct = np.rad2deg(phi_direct)
    theta = np.rad2deg(theta)
    
    obs3 = radiopropa.Observer()
    obsz3 = radiopropa.ObserverSurface(radiopropa.Sphere(radiopropa.Vector3d(channel[0], channel[1], channel[2]), 10 * radiopropa.meter))
    obs3.add(obsz3)
    obs3.setDeactivateOnDetection(True)
    sim.add(obs3)
    
    # Output
    output = radiopropa.TextOutput('output_traj.txt', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    #output.enable(radiopropa.Output.CurrentAmplitudeColumn)
    output.enable(radiopropa.Output.SerialNumberColumn)
    sim.add(output)
    
    sim.add(radiopropa.MaximumTrajectoryLength(1000*radiopropa.meter))
    
    #Start rays from 0 - 90 deg
    solutions_found = False
    number_of_solutions = 0

    for phi in np.arange(0,phi_direct, 1):
        x = hp.spherical_to_cartesian(np.deg2rad(x1_dir[0]), np.deg2rad(x1_dir[1]))
        y = hp.spherical_to_cartesian(np.deg2rad(phi), np.rad2deg(theta))
        delta = np.arccos(np.dot(x, y))
        
        cherenkov_angle = 56
        if (abs(np.rad2deg(delta) - cherenkov_angle) < 20): #only include rays with angle wrt cherenkov angle smaller than 20 degrees

            source = radiopropa.Source()
            
            source.add(radiopropa.SourcePosition(radiopropa.Vector3d(x1[0], x1[1], x1[2])))
            #source.add(radiopropa.SourceAmplitude(1))
            #source.add(radiopropa.SourceFrequency(1E6))
            #z = np.cos(phi * radiopropa.deg)
            #x = np.sin(phi * radiopropa.deg)
            x,y,z = hp.spherical_to_cartesian(phi * radiopropa.deg ,theta * radiopropa.deg)
            source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
            #print('Ray Direction {} deg = ({}, 0, {})'.format(phi, x, z))
            sim.setShowProgress(True)
            candidate = source.getCandidate()

            sim.run(candidate, True)
            trajectory_length = candidate.getTrajectoryLength()
            Candidate = candidate.get() #candidate is a pointer to the object Candidate
            
            detection = obsz3.checkDetection(Candidate) #return {0,1,2} voor {DETECTED,VETO,NOTHING}
            if detection == 0:
                pathx = np.fromstring(Candidate.getPathX()[1:-1],sep=',')
                pathy = np.fromstring(Candidate.getPathY()[1:-1],sep=',')
                pathz = np.fromstring(Candidate.getPathZ()[1:-1],sep=',')
                path = np.stack([pathx,pathy,pathz], axis=1)
                launchVector = Candidate.getLaunchVector()
                receiveVector = Candidate.getReceiveVector()
                print(launchVector,receiveVector)
