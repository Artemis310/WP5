import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import integrate
import matplotlib.pyplot as plt
import seaborn as sb

sb.set()
class AeroLoads:
    def __init__(self, rho, V_inf, alpha, dy):
        self.rho = rho
        self.V_inf = V_inf
        self.alpha = np.radians(alpha)
        self.dy = dy
        self.end = self.interpol()[-2][-1]
        self.num = int(self.end/self.dy)

    def interpol(self):
        la_0 = np.genfromtxt("MainWing_a=0.00_v=10.00ms.txt", skip_header=40, skip_footer=1029, encoding='iso-8859-1')
        la_10 = np.genfromtxt('MainWing_a=10.00_v=10.00ms.txt', skip_header=40, skip_footer=1029, encoding='iso-8859-1')

        ylst = []
        chrlst = []
        Cllst_a0 = []
        Cdlst_a0 = []
        Cmlst_a0 = []

        for i in range(len(la_0)):
            ylst.append(la_0[i][0])
            chrlst.append(la_0[i][1])
            Cllst_a0.append(la_0[i][3])
            Cdlst_a0.append(la_0[i][5])
            Cmlst_a0.append(la_0[i][7])

        Cllst_a10 = []
        Cdlst_a10 = []
        Cmlst_a10 = []

        for j in range(len(la_10)):
            Cllst_a10.append(la_10[j][3])
            Cdlst_a10.append(la_10[j][5])
            Cmlst_a10.append(la_10[j][7])

        chr_func = sp.interpolate.interp1d(ylst, chrlst, kind='linear', fill_value="extrapolate")
        Cl_func_a0 = sp.interpolate.interp1d(ylst, Cllst_a0, kind='cubic', fill_value="extrapolate")
        Cl_func_a10 = sp.interpolate.interp1d(ylst, Cllst_a10, kind='cubic', fill_value="extrapolate")
        Cd_func_a0 = sp.interpolate.interp1d(ylst, Cdlst_a0, kind='cubic', fill_value="extrapolate")
        Cd_func_a10 = sp.interpolate.interp1d(ylst, Cdlst_a10, kind='cubic', fill_value="extrapolate")
        Cm_func_a0 = sp.interpolate.interp1d(ylst, Cmlst_a0, kind='cubic', fill_value="extrapolate")
        Cm_func_a10 = sp.interpolate.interp1d(ylst, Cmlst_a10, kind='cubic', fill_value="extrapolate")

        return Cl_func_a0, Cl_func_a10, Cd_func_a0, Cd_func_a10, Cm_func_a0, Cm_func_a10, ylst, chr_func

    def dist(self, y):
        Cl_PerSpan_a0 = self.interpol()[0](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y)
        Cl_PerSpan_a10 = self.interpol()[1](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y)
        Cd_PerSpan_a0 = self.interpol()[2](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y)
        Cd_PerSpan_a10 = self.interpol()[3](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y)
        Cm_PerSpan_a0 = self.interpol()[4](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y) ** 2
        Cm_PerSpan_a10 = self.interpol()[5](y) * 1 / 2 * self.rho * self.V_inf ** 2 * self.interpol()[-1](y) ** 2

        return Cl_PerSpan_a0, Cl_PerSpan_a10, Cd_PerSpan_a0, Cd_PerSpan_a10, Cm_PerSpan_a0, Cm_PerSpan_a10

    def total_dist(self, load_fac = 2.5):
        param_sin = np.sin(self.alpha) / np.sin(np.radians(10))
        y = np.linspace(0, self.end, num=self.num)
        L_d = self.dist(y)[0] + param_sin * (self.dist(y)[1] - self.dist(y)[0])
        D_d = self.dist(y)[2] + param_sin * (self.dist(y)[3] - self.dist(y)[2])
        M_d = self.dist(y)[4] + param_sin * (self.dist(y)[5] - self.dist(y)[4])


        return y, load_fac * L_d, load_fac ** 2 * D_d, M_d, self.num

    def desired_dist(self, y):
        param_sin = np.sin(self.alpha) / np.sin(np.radians(10))
        L_d = self.dist(y)[0] + param_sin * (self.dist(y)[1] - self.dist(y)[0])
        D_d = self.dist(y)[2] + param_sin * (self.dist(y)[3] - self.dist(y)[2])
        M_d = self.dist(y)[4] + param_sin * (self.dist(y)[5] - self.dist(y)[4])

        return (y, "Span [m]"), (L_d, "Lift [N]"), (D_d, "Drag [N]"), (M_d, "Moment [Nm]")

    def plotting(self):
        plt.subplot(2, 2, 1)
        plt.title("L'(y)")
        plt.plot(self.total_dist()[0], self.total_dist()[1], 'b-')
        plt.xlabel("Span [m]")
        plt.ylabel("Lift per unit Span [N/m]")
        plt.grid(b = True, which = 'major')

        plt.subplot(2, 2, 2)
        plt.title("D'(y)")
        plt.plot(self.total_dist()[0], self.total_dist()[2], 'r-')
        plt.xlabel("Span [m]")
        plt.ylabel("Drag per unit Span [N/m]")
        plt.grid(b = True, which = 'major')

        plt.subplot(2, 2, 3)
        plt.title("M'(y)")
        plt.plot(self.total_dist()[0], self.total_dist()[3], 'g-')
        plt.xlabel("Span [m]")
        plt.ylabel("Moment per unit Span [N]")
        plt.grid(b = True, which = 'major')

        plt.tight_layout()
        plt.show()

    def find_total(self):
        L_tot_func = np.polyfit(self.total_dist()[0], self.total_dist()[1], 3)
        D_tot_func = np.polyfit(self.total_dist()[0], self.total_dist()[2], 3)
        M_tot_func = np.polyfit(self.total_dist()[0], self.total_dist()[3], 3)

        def L(y, param):
            return (L_tot_func[0] * y ** 3 + L_tot_func[1] * y ** 2 + L_tot_func[2] * y + L_tot_func[3]) * param

        def D(y, param):
            return (D_tot_func[0] * y ** 3 + D_tot_func[1] * y ** 2 + D_tot_func[2] * y + D_tot_func[3]) * param

        def M(y, param):
            return (M_tot_func[0] * y ** 3 + M_tot_func[1] * y ** 2 + M_tot_func[2] * y + M_tot_func[3]) * param

        L_tot = sp.integrate.quad(lambda y: L(y, 1), 0, self.end)[0]
        D_tot = sp.integrate.quad(lambda y: D(y, 1), 0, self.end)[0]
        M_tot = sp.integrate.quad(lambda y: M(y, 1), 0, self.end)[0]

        return L_tot, D_tot, M_tot, L, D, M


# (Density [kg/m3], Free Stream Velocity [m/s], Angle of Attack [degrees], dy]
# print(AeroLoads(0.333, 247.6, 4.42, 0.01).plotting())  # Uncomment for Plotting Total Distribution
# desired_dist([Input Span Location])
# print(AeroLoads(1.225, 10, 4.42, 0.01).desired_dist(10))  # Uncomment for Desired Value
# print(AeroLoads(0.333, 247.6, 4.42, 0.01).find_total()[0])  # Uncomment for Functions and Total Values
