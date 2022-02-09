import thermal as th
import matplotlib.pyplot as plt


Tw = []

Tw.append(th.equilibrium_wall_temperature_stagnation(40335, 3303, 0.034))

pos = [0.1, 0.25, 0.5, 1, 1.5, 3, 5]

for x in pos:
	Tw.append(th.equilibrium_wall_temperature_convective(40335, 3303, x, 5))

vec = [0, 0.1, 0.25, 0.5, 1, 1.5, 3, 5]
#print(Tw)

#plt.scatter(vec, Tw)

#plt.ylabel('Tw [K]')
#plt.xlabel('Axial Position [m]')
#plt.title('Radiative Equilibrium Tw, rn = 34 mm, 5 deg inclination')
#plt.grid('minor')
#plt.show()

print(th.double_integral_wavelength(7.8e-7, 1e-3, 0, 1))


