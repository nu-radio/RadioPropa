import radiopropa
import numpy as np



ice = radiopropa.ExponentialIndex(1.78, 0.51, 37.25)
firn_bottom = radiopropa.ExponentialIndex(1.68, 0.51, 37.25)
firn_top = radiopropa.ExponentialIndex(1.58, 0.51, 37.25)
iceModel = radiopropa.IceModel_DoubleFirn(firn_top, firn_bottom, ice, -20, -50, 0)

if __name__ == "__main__":
    # simulation setup
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, 1.))
    sim.add(radiopropa.MinimumAmplitude(.01))
    sim.add(radiopropa.MaximumTrajectoryLength(500))

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
    layer1 = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-100), radiopropa.Vector3d(0,0,1)),
                                      ice.getValue(radiopropa.Vector3d(0,0,-20)),
                                      firn_bottom.getValue(radiopropa.Vector3d(0,0,-20)))
    layer2 = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-50), radiopropa.Vector3d(0,0,1)),
                                      firn_bottom.getValue(radiopropa.Vector3d(0,0,-50)),
                                      firn_top.getValue(radiopropa.Vector3d(0,0,-50)))

    layer3 = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-30), radiopropa.Vector3d(0,0,1)),
                                      firn_bottom.getValue(radiopropa.Vector3d(0,0,-30)),
                                      firn_bottom.getValue(radiopropa.Vector3d(0,0,-30))-.2)
    layer4 = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-75), radiopropa.Vector3d(0,0,1)),
                                      ice.getValue(radiopropa.Vector3d(0,0,-75)),
                                      ice.getValue(radiopropa.Vector3d(0,0,-75))-.2)
    layer5 = radiopropa.Discontinuity(radiopropa.Plane(radiopropa.Vector3d(0,0,-15), radiopropa.Vector3d(0,0,1)),
                                      firn_top.getValue(radiopropa.Vector3d(0,0,-15)),
                                      firn_top.getValue(radiopropa.Vector3d(0,0,-15))-.2)
    sim.add(surface)
    sim.add(layer1)
    sim.add(layer2)
    sim.add(layer3)
    sim.add(layer4)
    sim.add(layer5)

    # Output
    output = radiopropa.HDF5Output('output_traj_45_bis.h5', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    #output.enable(radiopropa.Output.CurrentAmplitudeColumn)
    output.enable(radiopropa.Output.SerialNumberColumn)
    sim.add(output)

    # Source

    #Start rays from 0 - 90 deg
    for phi in [45]:

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

