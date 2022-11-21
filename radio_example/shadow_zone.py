import numpy as np
from copy import deepcopy as cp
import radiotools.helper as hp
from NuRadioReco.utilities import units
from NuRadioMC.utilities import medium
from matplotlib import pyplot as plt
import radiopropa as RP
import radiopropa as RP
import matplotlib.pyplot as plt
import numpy as np
import radiotools.helper as hp
from NuRadioReco.utilities import units
import h5py
from numpy import asarray
from numpy import savetxt
import csv
import pandas as pd


# Select data where x in (x_min, x_max) and z in (z_min, z_max)
x_min, x_max = 150, 150.4 	#Antenna Y position
z_min, z_max = -2, 0.1	#Antenna Z position
sourcepos = -37.5*RP.meter 	#Emitter Z position
reps = 10000 	#number of repetitions





#iceModel = medium.get_ice_model('greenland_simple').get_ice_model_radiopropa()
iceModel = medium.get_ice_model('ross_ice_shelf').get_ice_model_radiopropa()
#Air boundary for horizontal signals
airBoundary = RP.HorizontalSurface(RP.Plane(RP.Vector3d(0,0,0), RP.Vector3d(0,0,1)), 1.27, 1, surfacemode=False)


transmitter1 = np.array([0,0,-10*units.meter])
transmitter2 = np.array([0,0,-60*units.meter])
transmitter3 = np.array([0,0,-100*units.meter])



phi = 0*RP.deg		#scanning range minimum
step = 2*RP.deg		#scanning range steps
index=1
theta_scanning_range = np.arange(0,90*RP.deg+step, step)
for n in range(reps):
    for theta in theta_scanning_range: #below phi_direct no solutions are possible without upward reflections
        def radiopropa_simulation(trm_pos):
        ##define module list for simulation
            sim = RP.ModuleList()
            sim.add(RP.PropagationCK(iceModel.get_scalar_field(), 1E-8, .001, 1.))
            sim.add(airBoundary)


            #ADD SCATTERING LAYERS HERE: define sigma and layer depth and add layer to sim
            sigma = 0.05
            layerDepth = -15
            layerName = RP.RandomN(RP.Plane(RP.Vector3d(0,0,layerDepth), RP.Vector3d(0,0,1)),sigma)
            sim.add(layerName)



            #for module in iceModel.get_modules().values(): sim.add(module)
            sim.add(RP.MaximumTrajectoryLength(450*RP.meter))
            sim.add(RP.MinimumAmplitude(1e-9))

            ## define observer for stopping simulation (boundaries)
            obs = RP.Observer()
            obs.setDeactivateOnDetection(True)
            boundary_left = RP.ObserverSurface(RP.Plane(RP.Vector3d(450*RP.meter,0,0), RP.Vector3d(1,0,0)))
            obs.add(boundary_left)
            boundary_right = RP.ObserverSurface(RP.Plane(RP.Vector3d(-110*RP.meter,0,0), RP.Vector3d(1,0,0)))
            obs.add(boundary_right)
            boundary_bottom = RP.ObserverSurface(RP.Plane(RP.Vector3d(0,0,-110*RP.meter), RP.Vector3d(0,0,1)))
            obs.add(boundary_bottom)
            boundary_top = RP.ObserverSurface(RP.Plane(RP.Vector3d(0,0,5*RP.meter), RP.Vector3d(0,0,1)))
            obs.add(boundary_top)
            sim.add(obs)


            obs2 = RP.Observer()
            obs2.setDeactivateOnDetection(True)
            channel = RP.ObserverSurface(RP.Sphere(RP.Vector3d(400*RP.meter,0,-25*RP.meter),10*RP.meter))
            obs2.add(channel)
            sim.add(obs2)

            detected_rays = []

            ray = hp.spherical_to_cartesian(theta,phi)
            source = RP.Source()
            source.add(RP.SourcePosition(RP.Vector3d(0, 0,sourcepos)))
            source.add(RP.SourceAmplitude(1))
            source.add(RP.SourceFrequency(1E6))
            source.add(RP.SourceDirection(RP.Vector3d(*ray)))
            print('Ray Direction {} deg = ({}, 0, {})'.format(theta/RP.deg, ray[0], ray[2]))
            sim.setShowProgress(False)
            candidate = source.getCandidate()

            # Output
            output = RP.TextOutput('eventOutput.txt', RP.Output.Event3D)	#Output file name
            output.setLengthScale(RP.meter)
            output.enable(RP.Output.SerialNumberColumn)
            sim.add(output)


            sim.run(candidate, True)
        


        if __name__ == "__main__":
            radiopropa_simulation(transmitter3)

            '''
            rt = srt.ray_tracing(medium=medium.get_ice_model('greenland_simple'),shower_dir=-vertex_dir)
            rt.set_start_and_end_point(vertex_pos,channel_pos)
            rt.RadioPropa_raytracer()
            cand = rt._candidates
            rt.find_solutions()
            num = rt.get_number_of_solutions()
            
            print('len cand:',len(cand))
            print('num sol:',num)
            '''

	#FIND EVENTS:
        # Read txt file to pandas dataframe
        T_X_Z =  pd.read_csv('eventOutput.txt', delim_whitespace=True, comment="#", usecols=[1,5,7], names=["T","X","Z"])

        # Convert string NAN to np.nan, drop nan values
        T_X_Z = T_X_Z.replace(["-NAN","NAN"],np.nan).dropna()
        # Convert values from strings to floats, ignoring lines where errors occur
        T_X_Z = pd.DataFrame({"T": pd.to_numeric(T_X_Z["T"], errors = 'coerce'), "X": pd.to_numeric(T_X_Z["X"], errors = 'coerce'), "Z": pd.to_numeric(T_X_Z["Z"], errors = 'coerce')}).dropna()



        # Writing dataframe to file
        T_X_Z[["X","Z"]].to_csv("data_X_Y.csv",index=False, sep="\t", mode="a",header=False)

        #TXZ = np.array(T_X_Z)
        # Select data based on antenna position X & Z
        T_X_Z_sel = T_X_Z[(T_X_Z["X"] > x_min) & (T_X_Z["X"] < x_max) & (T_X_Z["Z"] > z_min) & (T_X_Z["Z"] < z_max)]


        # If at least one event was found (len(T_X_Z_sel != 0)), write time of that event to timedata.txt
        if len(T_X_Z_sel) != 0:

            T_X_Z_time = T_X_Z_sel["T"]
            #T_X_Z_sorted = T_X_Z_sel["T"].sort_values(ascending = True)
            # Scale time to ns
            end = len(T_X_Z_time)
            fileT1 = open("timedataStraight.txt","a")
            fileT2 = open("timedataHorizontal.txt", "a")
            
            if (end == 1):
                t_write = T_X_Z_time.min()*1e9
                fileT1.write(str(t_write) + "\n")
                print("Straight Event at "+ str(t_write) + " ns")

            #if theres more than one element, use the largest element and smallest element (air/firn pulse)
            if (end >1):
                #horizontal_sel
                t_write1 = T_X_Z_time.max()*1e9
                t_write2 = T_X_Z_time.min()*1e9
                #if (abs(t_write1 - t_write2) > 20) :
                fileT2.write(str(t_write1) + "\n")
                fileT2.write(str(t_write2) + "\n")
                print("Horizontal Event at "+ str(t_write2) + " ns and " + str(t_write1) + " ns" )
            # Open, write, close
            
            fileT1.close()
            fileT2.close()

            # Open, write, close
    print(index)
    index +=1
            
        #else: print("No event")


    