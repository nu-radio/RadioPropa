import radiopropa
import numpy as np


iceModel = radiopropa.GorhamIceModel(z0=-10.)
airBoundary = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,0), radiopropa.Vector3d(0,0,1)), 1.5 , 1.)
firnLayer = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-10.), radiopropa.Vector3d(0,0,1)), iceModel.getValue(radiopropa.Vector3d(0,0,-10.001)) , 1.5)


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
    for phi in np.linspace(0,90, 10):
        z = np.cos(phi * radiopropa.deg)
        x = np.sin(phi * radiopropa.deg)
        print('Ray Direction {} deg = ({}, 0, {})'.format(phi, x, z))
        source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
        sim.setShowProgress(True)
        sim.run(source, 1)
