import numpy as np
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad, dblquad
from scipy.optimize import fsolve
import constants as const
import skaero.atmosphere.coesa as sk

c = const.c # speed of light [m/s]
h = const.h # Planck's constant [m2 kg / s]
kb = const.kb # Boltzmann's constant [m2 kg s-2 K-1]

# aerothermal heating
# correlations due to Tauber, taken from Tracy & Wright

def stagnation_heating(Tw, y, v, rn):

	hw = v ** 2 + 2.3e5
	h0 = 1e3 * Tw
	
	q = 1.83e-4 / np.sqrt(rn) * (1 - (hw / h0)) * sk.density(y) * v ** 3

	return q

def convective_heating(Tw, y, v, pos, theta):

	hw = v ** 2 + 2.3e5
	h0 = 1e3 * Tw
	
	if v > 4e3:
		2.2e-5*(np.cos(np.radians(theta)) ** 2.08 * np.sin(np.radians(theta)) ** 1.6) / (pos ** 0.2) * (1 - 1.11 * hw / h0) * sk.density(y) ** 0.8 * v ** 3.7
	else:
		3.89e-4*(np.cos(np.radians(theta)) ** 1.78 * np.sin(np.radians(theta)) ** 1.6) / (pos ** 0.2) * (1 - 1.11 * hw / h0) * (556/Tw) ** 0.25 * sk.density(y) ** 0.8 * v ** 3.37
	
	return q

def equilibrium_wall_temperature():
	
	
	return Tw

# radiation emission due to surface heating

def spectral_radiance_wavelength(wavelength, Tw):
	""" Computes spectral radiance of emitted radiation as a function of wavelength and wall temperature.

	Args:
		wavelength: [m]
		wall temperature: [K]
	
	Returns:
		spectral radiance: [W m-3 sr-1]

	"""

	L = ((2 * h * c ** 2) / (wavelength ** 5)) * ( 1 / (np.exp(h * c / kb / Tw / wavelength) - 1))

	return L

def spectral_radiance_frequency(frequency, Tw):
	""" Computes spectral radiance of emitted radiation as a function of frequency and wall temperature.

	Args:
		frequency: [Hz]
		wall temperature: [K]
	
	Returns:
		spectral radiance: [W m-2 sr-1 Hz-1]

	"""

	L = ((2 * h * frequency ** 2) / (c ** 2)) * ( 1 / (np.exp(h * frequency / kb / Tw) - 1))

	return L

def temperature_distribution(x):
	""" TEST FUNCTION: Generates mock temperature distribution along 1D profile to test integration.

	Args:
		x: position relative to orign (positive) [m]
	
	Returns:
		T: temperature [K]

	"""

	T = 1000 - x*500

	return T

def double_integral_wavelength(wav_min, wav_max, x_min, x_max):
	""" Integrates over wavelength range and 1D region to find radiant intensity (radiant flux per solid angle element).

	Args:
		wav_min: wavelength lower limit on wavelength [m].
		wav_max: wavelength upper limit on wavelength [m].
	    x_min: lower limit on x [m].
		x_max: upper limit on x [m].

	Returns:
		radiant intensity [W sr-1]

	"""

    return dblquad(lambda x, wavelength: spectral_radiance_wavelength(wavelength, temperature_distribution(x)), wav_min, wav_max, x_min, x_max)

def double_integral_frequency(freq_min, freq_max, x_min, x_max):
	""" Integrates over frequency range and 1D region to find radiant intensity (radiant flux per solid angle element).

	Args:
		freq_min: wavelength lower limit on frequency [m].
		wav_max: wavelength upper limit on frequency [m].
	    x_min: lower limit on x [m].
		x_max: upper limit on x [m].

	Returns:
		radiant intensity [W sr-1]

	"""

    return dblquad(lambda x, wavelength: spectral_radiance_frequency(frequency, temperature_distribution(x)), freq_min, freq_max, x_min, x_max)

print(double_integral_wavelength(7.8e-7, 1e-3, 0, 1))
