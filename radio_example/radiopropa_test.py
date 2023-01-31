from NuRadioMC.SignalProp import propagation
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import NuRadioReco.framework.electric_field
import matplotlib.pyplot as plt
import numpy as np
import time


from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.SignalProp import analyticraytracing

pulse = np.array([[0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.,2.,0.,0.,0.,0.,0.]])
pulse = np.array([0.,0.,0.,0.,0.,1.,0.,0.,0.,0.])

delta = np.zeros(100000)
delta[50000] = 1


dt = 10**(-12) * units.second

sr = 1/dt

delta_pulse_ana = NuRadioReco.framework.electric_field.ElectricField([1], position=None,
                 shower_id=None, ray_tracing_id=None)
delta_pulse_num = NuRadioReco.framework.electric_field.ElectricField([1], position=None,
                 shower_id=None, ray_tracing_id=None)

delta_pulse_ana.set_trace(delta, sr)
delta_pulse_num.set_trace(delta, sr)

filter = delta_pulse_ana.get_filtered_trace([300* units.MHz, 700* units.MHz], filter_type = 'rectangular')
filter = 1/np.sqrt(2) * filter/max(filter)

zeros = np.zeros(100000)
delta_efield = np.vstack((filter, filter, zeros))


delta_pulse_ana.set_trace(delta_efield, sr)
delta_pulse_num.set_trace(delta_efield, sr)

solution_types = propagation.solution_types

# position a
initial_point = np.array([2500 , 0 , -1500 ])* units.m
final_point = np.array([0, 0, -200])* units.m

points = 100000


# ---------------analytical calculations---------------------
st = time.time()

ref_index_model = 'southpole_2015'
ref_index_model = 'birefringence_medium'
ice = medium.get_ice_model(ref_index_model)
rays = analyticraytracing.ray_tracing(ice)

rays.set_start_and_end_point(initial_point,final_point)
rays.find_solutions()
if 0:
    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[0], label='ana, theta, pre')
    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[1], '--', label='ana, phi, pre')
    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[2], ':', label='ana, R, pre')

    rays.apply_propagation_effects(delta_pulse_ana, 0)

    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[0], label='ana, theta, post')
    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[1], '--', label='ana, post')
    plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[2], ':', label='ana, post')

et = time.time()
elapsed_time = et - st
print('Execution time (analytical):', elapsed_time, 'seconds')


# ---------------radiopropa calculations---------------------
st = time.time()

ref_index_model = 'southpole_2015'
ice = medium.get_ice_model(ref_index_model)
rays = radioproparaytracing.radiopropa_ray_tracing(ice)

rays.set_start_and_end_point(initial_point,final_point)
rays.find_solutions()



if 0:
    print('test1')
    rays.apply_propagation_effects(delta_pulse_num, 0)
    print('test end')

if 1:
    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[0], label='num, theta, pre')
    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[1], '--', label='num, phi, pre')
    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[2], ':', label='num, R, pre')

    rays.apply_propagation_effects(delta_pulse_num, 0)

    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[0], label='num, theta, post')
    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[1], '--', label='num, post')
    plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace()[2], ':', label='num, post')

et = time.time()
elapsed_time = et - st

print('Execution time (numerical):', elapsed_time, 'seconds')

#plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[0], label='ana, theta')
#plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[1], '--', label='ana, phi')
#plt.plot(delta_pulse_ana.get_times(), delta_pulse_ana.get_trace()[2], ':', label='ana, R')


#plt.plot(delta_pulse_num.get_times(), delta_pulse_num.get_trace(), '--', label='num')
#plt.plot(ana_efield.get_trace())
#plt.plot(num_efield.get_trace())

#plt.xlim(450, 550)
plt.legend()
plt.show()




"""


for i_solution in range(rays.get_number_of_solutions()):

    solution_int = rays.get_solution_type(i_solution)
    solution_type = solution_types[solution_int]

    path = rays.get_path(i_solution, n_points= points)
    
    if i_solution == 0:
        path_a = path
    if i_solution == 1:
        path_b = path

end_a_x = path_a[:,0][-1]
end_a_y = path_a[:,2][-1]
end_b_x = path_b[:,0][-1]
end_b_y = path_b[:,2][-1]

dist = np.sqrt((end_a_x - end_b_x)**2 + (end_a_y - end_b_y)**2)     
print('endpoint distance: ' + str(round(dist, 4)) + 'm')


plt.scatter(path_a[:,0], path_a[:,2], label='refracted') 
plt.scatter(path_b[:,0], path_b[:,2], label='direct')
plt.legend(loc = 3) 
plt.show()

"""

