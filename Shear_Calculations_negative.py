import numpy as np
import scipy as sp
from scipy import integrate
import AeroLoads as ald #note that you might have to change this slightly
# depending on where you stored the Aeroloads file

import matplotlib.pyplot as plt

'''Constants'''
b_half = 51.73 / 2  # [m]
w_wing = 6566.81 * 9.81  # [N]
w_engine = 5942 * 9.81  # [N]
engine_location = 9.05  # [m]
engine_thrust = ((354 * 10 ** 3) / 2) * np.cos(np.radians(36.5))  # [N]
aerostandard = ald.AeroLoads(0.333, 247.66, 4.42, 0.001)

'''Distributed Loads'''
# Lift
L_fit = np.polyfit(aerostandard.total_dist()[0], 1.5*aerostandard.total_dist()[1], 4)

def lift(z):
    return -(-(L_fit[0] * z ** 4) - (L_fit[1] * z ** 3) - (L_fit[2] * z ** 2) - (L_fit[3] * z) - (
        L_fit[4]))                                                                                                     #Changed +-sign to --sign to account for -1g loading case.

def lift_centre_func(z):
    return ((L_fit[0] * z ** 4) + (L_fit[1] * z ** 3) + (L_fit[2] * z ** 2) + (L_fit[3] * z) + (
        L_fit[4])) * z

def centre_lift():
    return (sp.integrate.quad(lift_centre_func, 0, b_half)[0]) / (sp.integrate.quad(lift, 0, b_half)[0])

def LiftWeight(z):
    return -((L_fit[0] * z ** 4) + (L_fit[1] * z ** 3) + (L_fit[2] * z ** 2) + (L_fit[3] * z) + (
        L_fit[4])) - w_wing / b_half

# Wing weight
def DistributedWingWeight():
    return w_wing / b_half

# Drag
D_fit = np.polyfit(aerostandard.total_dist()[0], 1.5*aerostandard.total_dist()[2], 4)

def Drag(z):
    return (D_fit[0] * z ** 4) + (D_fit[1] * z ** 3) + (D_fit[2] * z ** 2) + (D_fit[3] * z) + (D_fit[4])

def drag_centre_func(z):
    return ((D_fit[0] * z ** 4) + (D_fit[1] * z ** 3) + (D_fit[2] * z ** 2) + (D_fit[3] * z) + (D_fit[4])) * z

def centre_drag():
    return (sp.integrate.quad(drag_centre_func, 0, b_half)[0]) / (sp.integrate.quad(Drag, 0, b_half)[0])


'''Point Loads'''
root_react_y = w_engine - sp.integrate.quad(LiftWeight, 0, b_half)[0]
root_react_x = sp.integrate.quad(Drag, 0, b_half)[0] - engine_thrust
dummy = sp.integrate.quad(LiftWeight, 0, b_half)[0]
bla = 9
'''Shear Distribution Function'''

def ShearDist(F_dist, F_root_react, F_engine, n_of_points):
    span = np.linspace(0, b_half, n_of_points)
    shear_value_at_span = np.zeros(shape=(len(span), 2))
    for i in range(n_of_points):
        shear_value_at_span[i, 0], shear_value_at_span[i, 1] = span[i], sp.integrate.quad(F_dist, 0, span[i])[
            0] + F_root_react + F_engine * np.heaviside(span[i] - engine_location, 0.5)
    return shear_value_at_span


'''Shear Distribution in yz plane'''
shear_yz = ShearDist(LiftWeight, root_react_y, -w_engine, 10000)
 # Plotting Shear Diagram
plt.plot(shear_yz[:, 0], shear_yz[:, 1])
plt.title("Shear Force Diagram for -1g")
plt.ylabel("Shear Force [N]")
plt.xlabel("Span [m]")
plt.show()

# Determining Shear Diagram Equation
engine_data_point_yz = np.where(shear_yz[:, 0] <= 9.05)[0][-1]
shear_yz_fit_before_engine = np.polyfit(shear_yz[:engine_data_point_yz + 1, 0], shear_yz[:engine_data_point_yz + 1, 1], 5)
shear_yz_fit_after_engine = np.polyfit(shear_yz[engine_data_point_yz:, 0], shear_yz[engine_data_point_yz:, 1], 5)

def ShearFuncBeforeEngineyz(z):
    return (shear_yz_fit_before_engine[0] * z ** 5) + (shear_yz_fit_before_engine[1] * z ** 4) + \
           (shear_yz_fit_before_engine[2] * z ** 3) + (shear_yz_fit_before_engine[3] * z ** 2) + \
           (shear_yz_fit_before_engine[4] * z) + (shear_yz_fit_before_engine[5])

def ShearFuncAfterEngineyz(z):
    return (shear_yz_fit_after_engine[0] * z ** 5) + (shear_yz_fit_after_engine[1] * z ** 4) + \
           (shear_yz_fit_after_engine[2] * z ** 3) + (shear_yz_fit_after_engine[3] * z ** 2) + \
           (shear_yz_fit_after_engine[4] * z) + (shear_yz_fit_after_engine[5])


'''Shear Distribution in xz plane'''
# Plotting Shear Diagram
shear_xz = ShearDist(lambda z: -Drag(z), root_react_x, engine_thrust, 10000)
# plt.plot(shear_xz[:, 0], shear_xz[:, 1])
# plt.ylabel("shear xz")
# plt.xlabel("span_test")
# plt.show()

# Determining Shear Diagram Equation
engine_data_point_xz = np.where(shear_xz[:, 0] <= 9.05)[0][-1]
shear_xz_fit_before_engine = np.polyfit(shear_xz[:engine_data_point_xz + 1, 0], shear_xz[:engine_data_point_xz + 1, 1], 5)
shear_xz_fit_after_engine = np.polyfit(shear_xz[engine_data_point_xz:, 0], shear_xz[engine_data_point_xz:, 1], 5)

def ShearFuncBeforeEnginexz(z):
    return (shear_xz_fit_before_engine[0] * z ** 5) + (shear_xz_fit_before_engine[1] * z ** 4) + \
           (shear_xz_fit_before_engine[2] * z ** 3) + (shear_xz_fit_before_engine[3] * z ** 2) + \
           (shear_xz_fit_before_engine[4] * z) + (shear_xz_fit_before_engine[5])

def ShearFuncAfterEnginexz(z):
    return (shear_xz_fit_after_engine[0] * z ** 5) + (shear_xz_fit_after_engine[1] * z ** 4) + \
           (shear_xz_fit_after_engine[2] * z ** 3) + (shear_xz_fit_after_engine[3] * z ** 2) + \
           (shear_xz_fit_after_engine[4] * z) + (shear_xz_fit_after_engine[5])