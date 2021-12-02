import scipy as sp
from scipy import integrate
import numpy as np
import AeroLoads as Al
import matplotlib.pyplot as plt
import matplotlib.axes
import TorsionalStiffness






t = 0.001  # [m] spar
s = 0.0016  # m flange
n = -1
s_factor = 1.5
thrust = np.cos(36.5) * (354 * (10 ** 3)) / 2  # [N]
engine_loc = 9.05  # [m]
engine_height = (4.979 - 0.89 - 1.524)  # [m]
G = 26 * 10 ** 9  # [Pa]
moment_func = Al.AeroLoads(0.333, 247.6, 4.42, 0.01).find_total()[5]
x_len = Al.AeroLoads(0.333, 247.6, 4.42, 0.01).total_dist()[0]
span = np.linspace(0, 51.73/2, 100)
span1 = np.linspace(0, engine_loc, 100)
span2 = np.linspace(engine_loc, 51.73/2 , 100)
num = Al.AeroLoads(0.333, 247.6, 4.42, 0.01).total_dist()[4]
engine_weight = 11884 * 9.80665 / 2  # [N]


def v_pos(x):
    y_pos = np.tan(np.radians(1.4)) * x + engine_height
    return y_pos
def h_pos(x):
    x_pos = 0.25 * (-0.2556 * x + 9.05)
    return x_pos
def polar(x, t, s):
    j = TorsionalStiffness.polar_moment_of_inertia(x, TorsionalStiffness.enclosed_area(x), t, s)
    return j

def torque_calc(y):
    moment_contr = sp.integrate.quad(lambda x: moment_func(x, n), 0, y)[0]
    thrust_contr = thrust * v_pos(engine_loc)
    engine_contr = engine_weight * h_pos(engine_loc)
    root_react_torque = engine_contr - thrust_contr - sp.integrate.quad(lambda x: moment_func(x, n), 0, 51.73/2)[0]
    total_torque = root_react_torque + moment_contr + (thrust_contr - engine_contr) * np.heaviside(y - engine_loc, 0.5)
    return total_torque * s_factor

torque_calc_vec = np.vectorize(torque_calc)
y1 = torque_calc_vec(span1)
y2 = torque_calc_vec(span2)

plt.plot(span, torque_calc_vec(span), label = "Torque [Nm]")
plt.title("Torque at spanwise location, n = -1")
plt.ylabel("Torque [Nm]")
plt.xlabel("Spanwise location[m]")
plt.ticklabel_format(axis = "y", style = "sci")
plt.legend()
plt.show()

polyfit_torque_b = np.polyfit(span1, y1, 4)
polyfit_torque_a = np.polyfit(span2, y2, 4 )

def torque_func_for_twist(y):
    return ((polyfit_torque_b[0] * y ** 4 + polyfit_torque_b[1] * y ** 3 + polyfit_torque_b[2] * y ** 2 + polyfit_torque_b[3] * y + polyfit_torque_b[4]) * np.heaviside(-y + engine_loc, 1) +
           (polyfit_torque_a[0] * y ** 4 + polyfit_torque_a[1] * y ** 3 + polyfit_torque_a[2] * y ** 2 + polyfit_torque_a[3] * y + polyfit_torque_a[4]) * np.heaviside(y - engine_loc, 1)) / (G * polar(y,t,s))

vec = np.vectorize(torque_func_for_twist)
# plt.plot(span, vec(span), label = "torque poly vs span")
# plt.legend()
# plt.show()

def find_twist(y):
        twist_dist = np.degrees(sp.integrate.quad(torque_func_for_twist , 0, y))[0]
        return twist_dist

find_twist_vec = np.vectorize(find_twist)

# plt.plot(span, find_twist_vec(span), label = "Twist")
# plt.title("Twist at spanwise location, n = -1")
# plt.xlabel("Spanwise location [m]")
# plt.ylabel("Twist [deg]")
# plt.legend()
# plt.show()


