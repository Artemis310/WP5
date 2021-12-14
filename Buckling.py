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


def corner_points(n_str_top, n_str_bot, width_str, 
                          area_str, centroid_x, centroid_y, th_spar, th_flang, height_str, thick, span_min = 0, span_max = 51.73/2):
    span_position = np.linspace(span_min, span_max, 100)
    
    NA_z = Mi.xx_vec_func(span_position, n_str_top, n_str_bot, width_str, 
                          area_str, centroid_x, centroid_y, th_spar, th_flang, height_str, thick)[0]
    Na_x = Mi.yy_vec_func(span_position, n_str_top, n_str_bot, width_str, 
                          area_str, centroid_x, centroid_y, th_spar, th_flang, height_str, thick)[0]

    top_right_z = Mi.y_coord1(Mi.c_spar2)*Mi.c_vec(span_position) - NA_z
    top_left_z = Mi.y_coord1(Mi.c_spar1)*Mi.c_vec(span_position) - NA_z
    bottom_z = (0.5*(Mi.y_coord2(Mi.c_spar1) + Mi.y_coord2(Mi.c_spar2)))*Mi.c_vec(span_position) - NA_z

    right_spar_x = Mi.c_spar2 * Mi.c(span_position) - Na_x
    left_spar_x = Mi.c_spar1 * Mi.c(span_position) - Na_x

    return top_right_z, top_left_z, bottom_z, right_spar_x, left_spar_x

corner_points_vec = np.vectorize(corner_points)


    
class NormalStressCalcs:
    def __init__(self, plane = None, cross_section_dist_z = 0, cross_section_dist_x = 0):
        self.plane = plane
        self.cross_section_dist_z = cross_section_dist_z
        self.cross_section_dist_x = cross_section_dist_x
        self.E = 69e9
        self.num = 100
    
    def stress_along_span(self, n_str_top, n_str_bot, width_str, 
                          area_str, centroid_x, centroid_y, th_spar, th_flang, height_str, thick, span_min = 0, span_max = 51.73/2):
        span_locations = np.linspace(span_min, span_max, self.num)

        inertia_xx = Mi.xx_vec_func(span_locations, n_str_top, n_str_bot, width_str, area_str, centroid_x, 
            centroid_y, th_spar, th_flang, height_str,thick)[1]
        inertia_yy = Mi.yy_vec_func(span_locations, n_str_top, n_str_bot, width_str, area_str, centroid_x, 
            centroid_y, th_spar, th_flang, height_str,thick)[1]

        if self.plane.lower == "Lift":
            return np.column_stack((span_locations, (Md.moment_yz_vec(span_locations) * self.cross_section_dist_z) / inertia_xx)) 
        elif self.plane.lower == "Drag":
            return np.column_stack((span_locations, (Md.moment_zx_vec(span_locations) * self.cross_section_dist_x) / inertia_yy))
        else:
            stress_values = (Md.moment_yz_vec(span_locations) * self.cross_section_dist_z) / inertia_xx
            + (Md.moment_zx_vec(span_locations) * self.cross_section_dist_x) / inertia_yy
            return np.column_stack((span_locations, stress_values))

    def find_stress_at_span(self, span_position):
        stress_index = np.where(self.stress_along_span()[:,0] <= span_position)[0][-1]
        return span_position, self.stress_along_span()[stress_index, 1], stress_index

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
        plt.ylabel("Normal Stress [Pa]")
        plt.grid(b = True, which = 'major')
        plt.legend()
        plt.show()

    def tension_analysis(self):
            plt.axhline(self.sigma_ult)
            self.plotting_stress()


class BuckleWeb:
    def __init__(self, t):
        self.t = t
        self.E = 69e9
        self.p_ratio = 0.33
        self.num = 1000
        self.span = np.linspace(0, 51.73 / 2, num=self.num)
        self.c_spar1 = Mi.c_spar1
        self.c_spar2 = Mi.c_spar2
        self.hf = np.linspace(Mi.t_c_spar1 * self.c_spar1 * Mi.c(0), Mi.t_c_spar1 * self.c_spar1 * Mi.c(self.span[-1]), num=self.num)
        self.hr = np.linspace(Mi.t_c_spar2 * self.c_spar2 * Mi.c(0), Mi.t_c_spar2 * self.c_spar2 * Mi.c(self.span[-1]), num=self.num)

    def cri_buckle_web(self, ks):
        tcr_f = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (self.t / self.hf) ** 2
        tcr_r = (np.pi**2 * ks * self.E) / (12 * (1 - self.p_ratio ** 2)) * (self.t / self.hr) ** 2

        return tcr_f, tcr_r

    def shear_ave(self):
        z = self.span
        V_yz = Sc.TotalShearyz(z)
        V_xz = Sc.TotalShearxz(z)
        shear_stress_yz = V_yz / (self.hf * self.t + self.hr * self.t)
        shear_stress_xz = V_xz / (self.hf * self.t + self.hr * self.t)

        return shear_stress_yz, shear_stress_xz

    def torque(self, z):
        return Tw.torque_calc_vec(z) / (2 * cross_section_area(z))  # Signs have to be checked

    def total_shear(self, ks):
        total = np.sqrt(self.shear_ave()[0]**2 + self.shear_ave()[1]**2) + self.torque(self.span)
        comparison = self.cri_buckle_web(ks)[0] - total

        if np.any(comparison <= 0):
            ans = False # Point(s) along the span have a higher stress than the critical
        else:
            ans = True # All is Good

        return total, ans, comparison

    def plotting_shear(self):
        plt.plot(self.span, self.total_shear(1)[0], 'r-' , label="yz-Plane")
        plt.plot(self.span, self.total_shear(1)[1], 'b-', label="xz-Plane")
        plt.xlabel("Span [m]")
        plt.ylabel("Shear Stress [Pa]")
        plt.grid(b = True, which = 'major')
        plt.legend()
        plt.show()

class Tension_analysis:
    def __init__(self, data_count= 1000, sigma_ult=310 * 10 ** 6, span_min= 0, span_max= 51.73 / 2):
        self.data_count = data_count
        self.span_locations = np.linspace(span_min, span_max, self.data_count)
        self.inertia_xx = Mi.xx_vec_func(self.span_locations, 2, 2, 0.1, 6 * 10 ** -4, 0.1, 0.1, 0.001, 0.001, 0.1, 0.002)
        self.inertia_yy = Mi.yy_vec_func(self.span_locations, 2, 2, 0.1, 6 * 10 ** -4, 0.1, 0.1, 0.001, 0.001, 0.1, 0.002)
        self.cross_section_dist_z_max = Mi.y_coord1(Mi.c_spar1) * Mi.c_vec(self.span_locations) - self.inertia_xx[0]
        self.cross_section_dist_x_max = Mi.c_spar1 * Mi.c_vec(self.span_locations) - self.inertia_yy[0]
        self.sigma_ult = sigma_ult

    def stress_along_span(self):
            second_tuple_val = (Md.moment_yz_vec(self.span_locations) * self.cross_section_dist_z_max) / \
                               self.inertia_xx[1] + (
                                           -Md.moment_zx_vec(self.span_locations) * self.cross_section_dist_x_max) / \
                               self.inertia_yy[1]
            return np.column_stack((self.span_locations, second_tuple_val)), second_tuple_val


    def plotting_stress(self):
        plot_label = "Stress due to yz and xz plane bending"
        plt.plot(self.stress_along_span()[:, 0], self.stress_along_span()[:, 1], 'k-', label=plot_label)
        plt.xlabel("Span [m]")
        plt.ylabel("Normal Stress [Pa]")
        plt.grid(b=True, which='major')
        plt.legend()
        plt.show()

    def tension_analysis(self):
        plt.axhline(self.sigma_ult)
        self.plotting_stress()




class BuckleSkin:
    def __init__(self, span_location, kc, t, stringer_count, stringer_width, plate_width):
        self.span_location = span_location
        self.kc = kc
        self.E = 69e9
        self.p_ratio = 0.33
        self.t = t
        self.stringer_count = stringer_count
        self.stringer_width = stringer_width
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

        max_stress_shear = max(BuckleWeb.total_shear[0], BuckleWeb.total_shear[1])

        fail_comp_normal = min(BuckleSkin(self.span_position).crit_buckle_skin, BuckleColumn(self.span_position).crit_buckle_stringer)
        fail_comp_shear = BuckleWeb(self.span_position).cri_buckle_web

        margin_of_safety_at_span = min(fail_comp_normal/max_stress_normal, fail_comp_shear/max_stress_shear)

        return margin_of_safety_at_span

    def plot_mos(self):
        plt.plot(self.span_position, self.find_mos[0])
        plt.xlabel("Span [m]")
        plt.ylabel("Margin of Safety [-]")
        plt.grid(b = True, which = 'major')
        plt.legend()
        plt.show()
        
class Design:
    def __init__(self, ks):
        self.ks = ks
    
    def buckle_check_web(self):
        i = 0.0001
        design_options = []
        while i < 0.05:
            if BuckleWeb(i).total_shear(self.ks)[1]:
                design_options.append(i)
                i +=  0.0001
            else:
                i +=  0.0001
        return design_options
    
    def buckle_check_skin(self):
        i = 0.0001
        BuckleSkin.crit_buckle_skin -Tension_analysis.stress_along_span()[1]


ks = 2
#print(Design(ks).buckle_check_web()[2])
