import radiopropa
import myPlugin

print("My Simulation\n")

ml = radiopropa.ModuleList()

ml.add(radiopropa.SimplePropagation(1*radiopropa.parsec, 100*radiopropa.parsec))
ml.add(radiopropa.MaximumTrajectoryLength(1000*radiopropa.parsec))
ml.add(myPlugin.MyModule())

print("+++ List of modules")
print(ml.getDescription())


print("+++ Preparing source")
source = radiopropa.Source()
source.add(myPlugin.AddMyProperty())
print(source.getDescription())

print("+++ Starting Simulation")
ml.run(source, 1)

print("+++ Done")
