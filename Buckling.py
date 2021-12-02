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

sns.set()

def stress_at_span(plane, span_location, cross_section_y):
    if plane.lower == "L":
        return (Md.moment_yz_vec(span_location)*cross_section_y) / Mi.moment_inertia_xx_func
    else:
        (Md.moment_zx_vec(span_location))*cross_section_y / Mi.moment_inertia_xx_func
    

def cross_section_area(y):
    return ((Mi.z1(y)+Mi.z4(y))*Mi.x3(y))/2

class BuckleWeb:
    def __init__(self, E, p_ratio,):
        self.E = E
        self.p_ratio = p_ratio
        self.span = np.linspace(0, 51.73 / 2, num=100)
        self.c_spar1 = Mi.c_spar1
        self.c_spar2 = Mi.c_spar2
        self.hf = np.linspace(Mi.t_c_spar1 * self.c_spar1 * Mi.c(0), Mi.t_c_spar1 * self.c_spar1 * Mi.c(self.span[-1]), num=100)
        self.hr = np.linspace(Mi.t_c_spar2 * self.c_spar2 * Mi.c(0), Mi.t_c_spar2 * self.c_spar2 * Mi.c(self.span[-1]), num=100)

    def cri_buckle_web(self, ks):
        tcr_f = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (0.1 / self.hf) ** 2
        tcr_r = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (0.1 / self.hr) ** 2

        return max(abs(tcr_f - tcr_r))

    def spar_geometry(self):
        tf = np.linspace(0.1, 0.5, num=100)
        tr = np.linspace(0.1, 0.5, num=100)

        return tf, tr

    def shear_ave(self):
        z = self.span
        V_yz = Sc.TotalShearyz(z)
        V_xz = Sc.TotalShearxz(z)
        shear_stress_yz = V_yz / (self.hf * self.spar_geometry()[0] + self.hr * self.spar_geometry()[1])
        shear_stress_xz = V_xz / (self.hf * self.spar_geometry()[0] + self.hr * self.spar_geometry()[1])

        return shear_stress_yz, shear_stress_xz

    def shear_flow(self):
        z = self.span
        T_yz = Sc.TotalShearyz(z)
        T_xz = Sc.TotalShearxz(z)
        A = cross_section_area(z)

        return T_yz / (2 * A), T_xz / (2 * A)

    def total_shear(self, ks):
        total_yz = (self.shear_ave()[0] + self.shear_flow()[0]) * 0.1
        comparison_yz = self.cri_buckle_web(ks) - total_yz
        total_xz = (self.shear_ave()[1] + self.shear_flow()[1]) * 0.1
        comparison_xz = self.cri_buckle_web(ks) - total_xz

        return total_yz, comparison_yz, total_xz, comparison_xz

    def plotting_shear(self):
        plt.plot(self.span, self.total_shear(1)[0], 'r-' , label="yz-Plane")
        plt.plot(self.span, self.total_shear(1)[2], 'b-', label="xz-Plane")
        plt.xlabel("Span [m]")
        plt.ylabel("Shear Stress [MPa]")
        plt.grid(b = True, which = 'major')
        plt.legend()
        plt.show()


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

print(BuckleWeb(69e9, 0.33).plotting_shear())
