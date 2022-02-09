import numpy as np
import atmosphere as atm
import constants as const
import matplotlib.pyplot as plt
import thermal as th

class dp_point:
    """ Class for holding extracted data points.

    """
    def __init__(self, x, y, v, rho, P, T, Ma, Tw):
        """ Initializer method.

        Args:
            x: range coordindate m
            y: altitude coordinate m
            v: velocity magnitude m/s
            rho: density kg/m3
            T: temperature K
            Ma: Mach number
            Tw: Equilibrium wall temperature K

        """
        self.x = x
        self.y = y
        self.v = v
        self.rho = rho
        self.P = P
        self.T = T
        self.Ma = Ma
        self.Tw = Tw

def find_nearest_index(val, target):
    """ Finds index of nearest solution points to queried values.

    Args:
        val: array-like list of values
        target: target value

    Returns:
        idx: index of closest value

    """

    array = np.asarray(val)
    idx = (np.abs(val - target)).argmin()
    
    return idx

def get_freestream_conditions(idx, fpa, x, y, v):
    """ Returns freestream conditions associated with index in trajectory solution.

    Args:
        idx: index to be queried
        fpa: array of flight path angles
        x: array of range 
        y: array of altitude
        v: array of velocity magnitude 

    Returns:
        dp_point object containing free-stream conditions

    """

    gamma = const.gamma
    R = const.R

    x_i = x
    y_i = y
    v_i = v
    rho_i = atm.get_density(y_i)
    P_i = atm.get_pressure(y_i)
    T_i = atm.get_temperature(y_i)
    Ma_i = v_i / np.sqrt(gamma * atm.get_R(y_i) * T_i)
    Tw_i = th.equilibrium_wall_temperature_convective(y, v, 1, 5)

    return dp_point(x_i, y_i, v_i, rho_i, P_i, T_i, Ma_i, Tw_i)

def process_solutions(dive_sol, glide_sol, vehicle, y0, v0):
    """ Converts solve_ivp solutions into arrays and concatenates arrays 
    from both dive and glide phases.

    Args:
        dive_sol: solution for dive phase
        glide_sol: solution for glide phase
        vehicle: vehicle_parameter object
        y0: altitude at start of glide phase
        v0: velocity at start of glide phase

    Returns:
        fpa: array of flight path angles
        x: array of range 
        y: array of altitude
        v: array of velocity magnitude 

    """

    L_D = vehicle.L_D
    b = vehicle.b
    g = const.g

    fpa = dive_sol.y[0]
    x = np.concatenate([dive_sol.y[1], glide_sol.y[0]])
    y = np.concatenate([dive_sol.y[2], glide_sol.y[1]])
    v = dive_sol.y[3]

    for i in range(len(glide_sol.y[0])):
        v_temp = np.sqrt(atm.get_density(y0)/atm.get_density(glide_sol.y[1][i])) * v0
        v = np.append(v, v_temp)
        fpa = np.append(fpa, 1 / vehicle.L_D / (1 + v_temp ** 2 / 2 / g / 8.4e3)) # TODO: variable scale height

    return fpa, x, y, v

def extract_trajectory_point_data(dive_sol, glide_sol, vehicle, y0, v0, pos_vec):
    """ Extra trajectory data at specified range positions.

    Args:
        dive_sol: solution for dive phase
        glide_sol: solution for glide phase
        vehicle: vehicle_parameter object
        y0: altitude at start of glide phase
        v0: velocity at start of glide phase
        pose_vec: vector of vehicle range positions

    Returns:
        points: list of point objects

    """

    fpa, x, y, v  = process_solutions(dive_sol, glide_sol, vehicle, y0, v0)

    points = []

    for pos in pos_vec:
        idx = find_nearest_index(x, pos)
        point = get_freestream_conditions(idx, fpa[idx], x[idx], y[idx], v[idx])
        points.append(point)

    return points
        

