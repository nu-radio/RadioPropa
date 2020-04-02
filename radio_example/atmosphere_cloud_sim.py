import radiopropa
import numpy as np


#Setup a simulation with source and observer in vacuum, layer with humid atmosphere above


#use implemented cloud model
cloud_h = 1000
#atmosphere as propagation medium with clouds from 1000 to 1500 m, weather at the ground: 300 Kelvin, 1000 mbar, 5% humidity
atmModel = radiopropa.CloudModel_atm(cloud_h, 1500, 300, 1000,0.05)
n_clouds = atmModel.getValue(radiopropa.Vector3d(0,0, cloud_h + 0.01))
n_air = atmModel.getValue(radiopropa.Vector3d(0,0, cloud_h - 0.01))
airBoundary = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0, cloud_h),radiopropa.Vector3d(0,0,1)), n_air, n_clouds) #air topped with clouds

if __name__=="__main__":
    #simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(atmModel, 1e-5, 0.1, 10))

    #create telescope as spherical observer
    obs = radiopropa.Observer()
    obs_surface = radiopropa.ObserverSurface(radiopropa.Sphere(radiopropa.Vector3d(17600,0,0),100))
    obs.add(obs_surface)
    obs.setDeactivateOnDetection(True)
    sim.add(obs)

    #cut off all rays going into earth, too high and behind detector
    end_planes = radiopropa.Observer()
    end_plane_bottom = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,0),radiopropa.Vector3d(0,0,-1)))
    end_plane_top = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,10000),radiopropa.Vector3d(0,0,1)))
    end_plane_right = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(20000,0,0),radiopropa.Vector3d(1,0,0)))
    end_planes.add(end_plane_bottom)
    end_planes.add(end_plane_top)
    end_planes.add(end_plane_right)
    end_planes.setDeactivateOnDetection(True)
    sim.add(end_planes)

    #output
    output = radiopropa.HDF5Output('trajectory_atmosphere.h5', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    output.enable(radiopropa.Output.SerialNumberColumn)
    sim.add(output)

    #Simulate rays from 0 to 90 degrees in 5 degree steps
    angles = np.arange(0,90,10)
    for phi in angles:
        sim.add(airBoundary)
        sim.add(radiopropa.MaximumTrajectoryLength(50000*radiopropa.meter))

        source = radiopropa.Source()
        source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0,0,0)))
        source.add(radiopropa.SourceAmplitude(1))
        source.add(radiopropa.SourceFrequency(50e6))

        z = np.cos(phi * radiopropa.deg)
        x = np.sin(phi * radiopropa.deg)
        source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
        print('Ray Direction {} deg').format(phi)
        sim.setShowProgress(True)
        sim.run(source, 1)


