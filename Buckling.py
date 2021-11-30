import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import integrate
from scipy import interpolate
from WP_4 import AeroLoads as Al
from WP_4 import Shear_Calculations as Sc


class BuckleWeb:
    def init(self, E, p_ratio, t, web_width, web_height, spar_height_front, spar_t_front, spar_height_rear,
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
        Vx = Sc.TotalShearDistxz(z)
        Vy = Sc.TotalShearDistyz(z)
        shear_stress_x = Vx / (self.hf * self.tf + self.hr * self.tr)
        shear_stress_y = Vy / (self.hf * self.tf + self.hr * self.tr)

        return shear_stress_x, shear_stress_y
class BuckleSkin:
    def __init__(self, kc ,E, t, stringer_count, stringer_width, p_ratio):
        self.kc = kc
        self.E = E
        self.t = t
        self.stringer_count = stringer_count
        self.stringer_width = stringer_width
        self.p_ratio = p_ratio

    def crit_buckle_skin(self):
        return  (((np.pi**2)*self.kc*self.E) / (12 * (1-self.p_ratio**2))) * ((self.t / self.t))**2