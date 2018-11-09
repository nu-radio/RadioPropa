import radiopropa

n0 = 2.
a = 1.
iceModel = radiopropa.n2linear(n0, a)

if __name__ == "__main__":
    sim = radiopropa.ModuleList()
    sim.add(radiopropa.PropagationCK(iceModel, 1E-8, .001, .01))

    output = radiopropa.HDF5Output('output_traj.h5', radiopropa.Output.Trajectory3D)
    output.setLengthScale(radiopropa.meter)
    sim.add(output)

    source = radiopropa.Source()
    source.add(radiopropa.SourcePosition(radiopropa.Vector3d(0, 0, 0.)))
    source.add(radiopropa.SourceParticleType(radiopropa.nucleusId(1, 1)))
    source.add(radiopropa.SourceAmplitude(1))
    source.add(radiopropa.SourceFrequency(1E6))


    sim.add(radiopropa.MaximumTrajectoryLength(100 * radiopropa.meter))

    sim.setShowProgress(True)
    source.add(radiopropa.SourceDirection(radiopropa.Vector3d(1., 0 , 1.)))
    sim.run(source, 1)
