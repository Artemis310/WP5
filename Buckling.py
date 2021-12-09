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
import TwistDistribution_requires_editing as Tw

sns.set()
#Material Constants
tens_yield_strength = 276e6
ult_yield_strength = 310e6

def cross_section_area(y):
    return ((Mi.z1(y)+Mi.z4(y))*Mi.x3(y))/2


def corner_points(span_position):
    NA_z = Mi.moment_inertia_xx_func()[0]
    Na_x = Mi.moment_inertia_yy_func()[0]

    top_right_z = Mi.y_coord1(span_position)*0.61*Mi.c(span_position) - NA_z
    top_left_z = Mi.y_coord1(span_position)*0.2*Mi.c(span_position) - NA_z
    bottom_z = Mi.y_coord2(span_position)*0.61*Mi.c(span_position) - NA_z

    right_spar_x = 0.61*Mi.c(span_position) - Na_x
    left_spar_x = 0.2*Mi.c(span_position) - Na_x

    return top_right_z, top_left_z, bottom_z, right_spar_x, left_spar_x

corner_points_vec = np.vectorize(corner_points)

    
class NormalStressCalcs:
    def __init__(self, plane = None, cross_section_dist_z = 0, cross_section_dist_x = 0):
        self.plane = plane
        self.cross_section_dist_z = cross_section_dist_z
        self.cross_section_dist_x = cross_section_dist_x
        self.E = 69e9
        self.num = 100
    
    def stress_along_span(self, span_min = 0, span_max = 51.73/2):
        span_locations = np.linspace(span_min, span_max, self.num)

        if self.plane.lower == "Lift":
            return np.column_stack((span_locations,(Md.moment_yz_vec(span_locations)*self.cross_section_dist_z) / 1)) #Mi.moment_inertia_xx_func(span_location)))
        elif self.plane.lower == "Drag":
            return np.column_stack((span_locations, (Md.moment_zx_vec(span_locations)*self.cross_section_dist_x) / 1)) #Mi.moment_inertia_yy_func(span_location)))
        else:
            return np.column_stack((span_locations, (Md.moment_yz_vec(span_locations)*self.cross_section_dist_z) / 1
            + (Md.moment_zx_vec(span_locations)*self.cross_section_dist_x) / 1))  #double check this equation, as well as the one above

    def find_stress_at_span(self, span_position):
        stress_index = np.where(self.stress_along_span()[:,0] <= span_position)[0][-1]
        return span_position, self.stress_along_span()[stress_index, 1]

    def compres_failure(self):
        yield_stress = 276e6 # Might be incorrect still has to be checked later
        total_len = np.linspace(0, 51.73/2, self.num)
        cond_1 = NormalStressCalcs("Combined", corner_points(total_len)[0], corner_points(total_len)[3]).find_stress_at_span(total_len)
        cond_2 = NormalStressCalcs("Combined", corner_points(total_len)[1], corner_points(total_len)[-1]).find_stress_at_span(total_len)
        cond_3 = NormalStressCalcs("Combined", corner_points(total_len)[2], corner_points(total_len)[-1]).find_stress_at_span(total_len)
        cond_4 = NormalStressCalcs("Combined", corner_points(total_len)[2], corner_points(total_len)[3]).find_stress_at_span(total_len)
        if cond_1.any() or cond_2.any() or cond_3.any() or cond_4.any() >= yield_stress:
            ans = "Point(s) along the span have a higher stress than the yield stress"
        else:
            ans = "All is Good"
        return ans

    def plotting_stress(self):
        if self.plane.lower == "Lift":
            plot_label = "yz Plane"
        elif self.plane.lower == "Drag":
            plot_label = "xz Plane"
        else:
            plot_label = "Stress due to yz and xz plane bending"
        plt.plot(self.stress_along_span()[:,0], self.stress_along_span()[:,1], 'k-' , label = plot_label)
        plt.xlabel("Span [m]")
        plt.ylabel("Normal Stress [MPa]")
        plt.grid(b = True, which = 'major')
        plt.legend()
        plt.show()


class BuckleWeb:
    def __init__(self):
        self.E = 69e9
        self.p_ratio = 0.33
        self.num = 100
        self.span = np.linspace(0, 51.73 / 2, num=self.num)
        self.c_spar1 = Mi.c_spar1
        self.c_spar2 = Mi.c_spar2
        self.hf = np.linspace(Mi.t_c_spar1 * self.c_spar1 * Mi.c(0), Mi.t_c_spar1 * self.c_spar1 * Mi.c(self.span[-1]), num=self.num)
        self.hr = np.linspace(Mi.t_c_spar2 * self.c_spar2 * Mi.c(0), Mi.t_c_spar2 * self.c_spar2 * Mi.c(self.span[-1]), num=self.num)

    def spar_geometry(self):
        tf = np.linspace(0.1, 0.5, num=self.num)
        tr = np.linspace(0.1, 0.5, num=self.num)

        return tf, tr

    def cri_buckle_web(self, ks):
        tcr_f = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (self.spar_geometry()[0] / self.hf) ** 2
        tcr_r = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (self.spar_geometry()[1] / self.hr) ** 2

        return tcr_f, tcr_r

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

    def torque(self, z):
        return Tw.torque_calc_vec(z) / (2 * cross_section_area(z))  # Signs have to be checked

    def total_shear(self, ks=2):
        total_yz = (self.shear_ave()[0] + self.shear_flow()[0] + self.torque(self.span)) * 0.1 
        comparison_yz = self.cri_buckle_web(ks)[0] - total_yz
        total_xz = (self.shear_ave()[1] + self.shear_flow()[1] + self.torque(self.span)) * 0.1
        comparison_xz = self.cri_buckle_web(ks)[1] - total_xz

        if comparison_yz.any() or comparison_xz.any() > 0:
            ans = "Point(s) along the span have a higher stress than the critical"
        else:
            ans = "All is Good"

        return total_yz, total_xz, ans

    def plotting_shear(self):
        plt.plot(self.span, self.total_shear(1)[0], 'r-' , label="yz-Plane")
        plt.plot(self.span, self.total_shear(1)[1], 'b-', label="xz-Plane")
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
class MarginOfSafety:
    def __init__(self, span_position):
        self.span_position = span_position

    def find_mos(self):
        applied_stress_top_right = NormalStressCalcs("Combined", corner_points_vec(self.span_position)[0], corner_points_vec(self.span_position)[3]).find_stress_at_span(self.span_position)
        applied_stress_top_left = NormalStressCalcs("Combined", corner_points_vec(self.span_position)[1], corner_points_vec(self.span_position)[-1]).find_stress_at_span(self.span_position)
        applied_stress_bottom_left = NormalStressCalcs("Combined", corner_points_vec(self.span_position)[2], corner_points_vec(self.span_position)[-1]).find_stress_at_span(self.span_position)
        applied_stress_bottom_right = NormalStressCalcs("Combined", corner_points_vec(self.span_position)[2], corner_points_vec(self.span_position)[3]).find_stress_at_span(self.span_position)

        max_stress_normal = max(applied_stress_top_right, applied_stress_top_left, applied_stress_bottom_right, applied_stress_bottom_left)

        fail_comp_normal = min(BuckleSkin(self.span_position).crit_buckle_skin, BuckleColumn(self.span_position).crit_buckle_stringer)
        fail_comp_shear = BuckleWeb(self.span_position).cri_buckle_web

        margin_of_safety_at_span = min(fail_comp_normal/max_stress_normal, fail_comp_shear/BuckleWeb.total_shear[0])

        return margin_of_safety_at_span

    def plot_mos(self):
        
        return None
        

<<<<<<< HEAD
print(BuckleWeb().plotting_shear())
=======

ks = 2
print(BuckleWeb().plotting_shear(), BuckleWeb().total_shear(ks)[2])
>>>>>>> 86a09df3207131d1d0f99a323d23be26885c8721


#print(StressCalcs("lift", 0.5, 1000).stress_along_span())
# stress = StressCalcs("Lift", 0.5, 0, 1000)

<<<<<<< HEAD
# stress.plotting_stress()
=======
# stress.plotting_stress()
>>>>>>> 86a09df3207131d1d0f99a323d23be26885c8721
