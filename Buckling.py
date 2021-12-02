import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import integrate
from scipy import interpolate
from WP_4 import AeroLoads as Al
from WP_4 import Shear_Calculations as Sc
from WP_4 import moment_inerta as Mi
from WP_4 import moment_diagram



def stress_at_span(span_pos, cross_section_y):
    return 
    

def cross_section_area(y):
    return ((Mi.z1+Mi.z4)*Mi.x3)/2

class BuckleWeb:
    def __init__(self, E, p_ratio, t, web_width, web_height, spar_height_front, spar_t_front, spar_height_rear,
                 spar_t_rear):
        self.E = E
        self.p_ratio = p_ratio
        self.t = t
        self.b = web_width
        self.h = web_height
        self.hf = spar_height_front
        self.tf = spar_t_front
        self.hr = spar_height_rear
        self.tr = spar_t_rear
        self.span = Al.AeroLoads(0.333, 247.66, 4.42, 0.01).total_dist()[0]

    def cri_buckle_web(self, ks):
        return (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio**2)) * (self.t / self.b) ** 2

    def shear_ave(self):
        z = self.span
        V = Sc.TotalShearDistxz(z)
        shear_stres = V / (self.hf * self.tf + self.hr * self.tr)

        return shear_stres

    def shear_flow(self):
        z = self.span
        T = Sc.TotalShearDistxz(z)
        A = cross_section_area(z)

        return T / (2 * A)

    def total_shear(self):
        return self.shear_ave + self.shear_flow * self.t

    def plotting_shear(self):
        plt.plot(self.span, self.total_shear)
        plt.xlabel("Span [m]")
        plt.ylabel("Shear Stress [MPa]")
        plt.plot()


class BuckleSkin:
    def __init__(self, span_location, kc ,E, t, stringer_count, stringer_width, p_ratio, plate_width):
        self.span_location = span_location
        self.kc = kc
        self.E = E
        self.t = t
        self.stringer_count = stringer_count
        self.stringer_width = stringer_width
        self.p_ratio = p_ratio
        self.plate_width = plate_width
        self.b = plate_width/self.stringer_count - stringer_width

    def crit_buckle_skin(self):
        return  (((np.pi**2)*self.kc*self.E) / (12 * (1-self.p_ratio**2))) * ((self.t / self.t))**2

class BuckleColumn:
    def __init__(self, K, E, I, L, A):
        self.K = K
        self.E = E
        self.I = I
        self.L = L
        self.A = A

    def crit_buckle_stringer(self):
        return (self.K * np.pi**2 * self.E * self.I) / (self.L**2 * self.A)
