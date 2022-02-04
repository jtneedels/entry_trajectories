import skaero.atmosphere.coesa as sk
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import constants as const
from dynamics import dive_pull_dynamics, glide_dynamics
from vehicle_parameters import vehicle_parameters
from data_processing import process_solutions, extract_trajectory_point_data

# define vehicle
missile = vehicle_parameters(2.562, 2992, 0, 0)

# entry parameters
v_o = 5e3 # m/s
h_o = 80e3 # m
fpa_o_vec = [np.radians(5)] # radians
q_lim_vec = [20e3] # Pa
alt_lim = 100 # m

for fpa_o in fpa_o_vec:
	for q_lim in q_lim_vec:

		print(fpa_o, q_lim)

		def dynamic_pressure_limit(t, z, vehicle):

			fpa, x, y, v = z

			return .5 * sk.density(y) * v ** 2 - q_lim

		def altitude_limit(t, z, v0, y0, vehicle):

			x, y = z

			return y - 100

		dynamic_pressure_limit.terminal = True 
		dynamic_pressure_limit.direction = -1 

		altitude_limit.terminal = True 

		dive_sol = solve_ivp(dive_pull_dynamics, [0, 5e3], [fpa_o, 0, h_o, v_o], method='RK45', events=dynamic_pressure_limit, args=[missile], max_step=1)

		y0 = dive_sol.y_events[0][0][2]
		v0 = dive_sol.y_events[0][0][3]

		glide_sol = solve_ivp(glide_dynamics, [0, 5e3], [dive_sol.y_events[0][0][1], dive_sol.y_events[0][0][2]], method='RK45', events=altitude_limit, args=(v0, y0, missile), max_step=1)

		plt.plot(np.concatenate([dive_sol.y[1]/1e3, glide_sol.y[0]/1e3]), np.concatenate([dive_sol.y[2]/1e3, glide_sol.y[1]/1e3]), label='V0: '+str(int(v_o/1e3))+'kps, EFPA: '+str(int(np.degrees(fpa_o)))+'deg, q: '+str(int(q_lim/1e3))+'kPa')

plt.xlabel('Distance [km]')
plt.ylabel('Altitude [km]')
plt.title('Trajectory')
plt.legend()
plt.grid('minor')

fpa, x, y, v = process_solutions(dive_sol, glide_sol, missile, y0, v0)
points = extract_trajectory_point_data(dive_sol, glide_sol, missile, y0, v0, [100e3, 200e3, 300e3, 400e3, 500e3, 1000e3, 1500e3, 2000e3, 2500e3])

for point in points:
	print('Alt: ', point.y, 'Press: ', point.P, 'Temp: ', point.T, 'Ma: ', point.Ma)
	plt.scatter(point.x/1e3, point.y/1e3, color='black')

plt.show()
