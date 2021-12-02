import numpy as np
import matplotlib.pyplot as plt
import Shear_Calculations as scalc
import scipy as sp
from scipy import integrate


# Inputs
Engine_thrust = scalc.engine_thrust
Engine_thrust_loc = scalc.engine_location
Engine_weight = scalc.w_engine
wing_weight = scalc.w_wing
Engine_loc = scalc.engine_location
Lift = sp.integrate.quad(scalc.lift, 0, scalc.b_half)[0]
Drag = sp.integrate.quad(scalc.Drag, 0, scalc.b_half)[0]
loc_lift = scalc.centre_lift()
loc_drag = scalc.centre_drag()
quarter_span = scalc.b_half / 2
span = np.linspace(0, 25.865, 1000)

# Equilibrium Equations
M1 = - Engine_weight * Engine_loc + Lift * loc_lift - wing_weight * quarter_span
M2 = - Drag * loc_drag + Engine_thrust * Engine_loc

def Moment_diagramzy(y):
    Mom_Func_ZY = M1 + (sp.integrate.quad(scalc.ShearFuncBeforeEngineyz, 0, y)[0]) * np.heaviside((-y + Engine_loc), 0) + \
                   (sp.integrate.quad(scalc.ShearFuncBeforeEngineyz, 0, Engine_loc)[0]) * np.heaviside((y - Engine_loc), 0) + \
                  (sp.integrate.quad(scalc.ShearFuncAfterEngineyz, Engine_loc, y)[0]) * np.heaviside((y - Engine_loc), 0)
    return Mom_Func_ZY

def moment_diagram_zx(y):
    mom_func_zx = M2 + (sp.integrate.quad(scalc.ShearFuncBeforeEnginexz, 0, y)[0]) * np.heaviside((-y + Engine_loc), 0) + \
                  (sp.integrate.quad(scalc.ShearFuncBeforeEnginexz, 0, Engine_loc)[0]) * np.heaviside((y - Engine_loc),0) + \
                  (sp.integrate.quad(scalc.ShearFuncAfterEnginexz, Engine_loc, y)[0]) * np.heaviside((y - Engine_loc),0)
    return mom_func_zx

# #vectorization of functions
moment_yz_vec = np.vectorize(Moment_diagramzy)
moment_zx_vec = np.vectorize(moment_diagram_zx)

""" Uncomment here under to get a plot of the moment diagrams of both planes"""
# #computation and plotting
# moment_yz = moment_yz_vec(span)
# moment_zx = moment_zx_vec(span)
# plt.plot(span, moment_yz,label = "yz")
# plt.plot(span, moment_zx, label = "zx")
# plt.ylabel("Moment zx-yz")
# plt.xlabel("span")
# plt.legend()
# plt.show()
