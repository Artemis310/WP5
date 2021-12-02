import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import integrate
from scipy import interpolate
import AeroLoads as Al
import Shear_Calculations as Sc
import moment_inerta as Mi
import moment_diagram as Md



def cross_section_area(y):
    return ((Mi.z1+Mi.z4)*Mi.x3)/2

class StressCalcs:
    def __init__(self, plane, cross_section_y, data_count):
        self.plane = plane
        self.cross_section_y = cross_section_y
        self.data_count =  data_count
    
    def stress_along_span(self):
        span_location = np.linspace(0, 51.73/2, self.data_count)
    
        if self.plane.lower == "Lift":
            return np.column_stack((span_location,(Md.moment_yz_vec(span_location)*self.cross_section_y) / Mi.moment_inertia_xx_func(span_location)))
        else:
            return np.column_stack((span_location, (Md.moment_zx_vec(span_location)*self.cross_section_y) / Mi.moment_inertia_xx_func(span_location)))

    def find_stress_at_span(self, span_position):
        return np.where(self.stress_along_span()[:,1] <= span_position)

    def plotting_stress(self):
        return None


class BuckleWeb:
    def __init__(self):
        self.E = 69e3
        self.p_ratio = 0.33
        self.span = np.linspace(0, 51.73 / 2, num=100)
        self.c_spar1 = Mi.c_spar1
        self.c_spar2 = Mi.c_spar2
        self.hf = np.linspace(Mi.t_c_spar1 * self.c_spar1 * Mi.c(0), Mi.t_c_spar1 * self.c_spar1 * Mi.c(self.span[-1]), num=100)
        self.hr = np.linspace(Mi.t_c_spar2 * self.c_spar2 * Mi.c(0), Mi.t_c_spar2 * self.c_spar2 * Mi.c(self.span[-1]), num=100)

    def cri_buckle_web(self, ks):
        tcr_f = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (0.1 / self.hf) ** 2
        tcr_r = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (0.1 / self.hr) ** 2

        return tcr_f, tcr_f

    def spar_geometry(self):
        tf = np.linspace(0.1, 0.5, num=100)
        tr = np.linspace(0.1, 0.5, num=100)

        return tf, tr

    def shear_ave(self):
        z = self.span
        V = Sc.TotalShearDistxz(z)
        shear_stres = V / (self.hf * self.tf + self.hr * self.tr)

        return shear_stres

    def shear_flow(self):
        z = self.span
        T = Sc.TotalShearDistxz(z)
        A = cross_section_area(z)

        return T_yz / (2 * A), T_xz / (2 * A)

    def total_shear(self, ks):
        total_yz = (self.shear_ave()[0] + self.shear_flow()[0]) * 0.1
        comparison_yz = self.cri_buckle_web(ks)[0] - total_yz
        total_xz = (self.shear_ave()[1] + self.shear_flow()[1]) * 0.1
        comparison_xz = self.cri_buckle_web(ks)[1] - total_xz

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

print(BuckleWeb().plotting_shear(), BuckleWeb().total_shear()[1:3:1])
