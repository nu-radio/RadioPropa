import radiopropa
import numpy as np



ice = radiopropa.ExponentialIndex(1.78, 0.51, 37.25)
firn_bottom = radiopropa.ExponentialIndex(1.58, 0.51, 37.25)
firn_top = radiopropa.ExponentialIndex(1.38, 0.51, 37.25)
iceModel = radiopropa.IceModel_DoubleFirn(firn_top, firn_bottom, ice, -20, -50, 0)

if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))

    # Observer to stop imulation at =0m and z=300m
    obs = radiopropa.Observer()
    obsz = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,10), radiopropa.Vector3d(0,0,1)))
    obs.add(obsz)
    obs.setDeactivateOnDetection(True)
    sim.add(obs)

    obs2 = radiopropa.Observer()
    obsz2 = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0,0,-300), radiopropa.Vector3d(0,0,1)))
    obs.add(obsz2)
    obs2.setDeactivateOnDetection(True)
    sim.add(obs2)

    surface = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,0), radiopropa.Vector3d(0,0,1)),
                                       firn_top.getValue(radiopropa.Vector3d(0,0,0)),1.)
    firn_top = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-20), radiopropa.Vector3d(0,0,1)),
                                        firn_bottom.getValue(radiopropa.Vector3d(0,0,-20)),
                                        firn_top.getValue(radiopropa.Vector3d(0,0,-20)))
    firn_bottom = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-50), radiopropa.Vector3d(0,0,1)),
                                           ice.getValue(radiopropa.Vector3d(0,0,-50)),
                                           firn_bottom.getValue(radiopropa.Vector3d(0,0,-50)))
    sim.add(surface)
    #sim.add(firn_top)
    #sim.add(firn_bottom)

    # Output
    output = radiopropa.HDF5Output('output_traj.h5', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    #output.enable(radiopropa.Output.CurrentAmplitudeColumn)
    output.enable(radiopropa.Output.SerialNumberColumn)
    #sim.add(output)

    # Source

    #Start rays from 0 - 90 deg
    for phi in [30]:

        source = radiopropa.Source()
        source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0, 0, -240.)))
        source.add(radiopropa.SourceAmplitude(1))
        source.add(radiopropa.SourceFrequency(1E6))
        z = np.cos(phi * radiopropa.deg)
        x = np.sin(phi * radiopropa.deg)
        source.add(radiopropa.SourceDirection(radiopropa.Vector3d(x, 0 , z)))
        print('Ray Direction {} deg = ({}, 0, {})'.format(phi, x, z))
        sim.setShowProgress(True)
        sim.run(source, 1)

    