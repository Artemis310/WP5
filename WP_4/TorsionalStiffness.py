from moment_inerta import z1,z4,x1,x3, y_coord1, c_spar1, c_spar2, c
import matplotlib.pyplot as plt
G = 26*10**9

def enclosed_area(y):
    return z4(y) * x3(y) + 1/2 * x3(y) * (y_coord1(c_spar1) - y_coord1(c_spar2)) * c(y)
def polar_moment_of_inertia(y, S, th_spar, th_flang):
    J = 4*S**2 / ((x3(y) + x1(y))/th_flang + (z1(y) + z4(y))/th_spar)
    return J
def torsional_stiffness(J,G):
    t_s = J*G
    return t_s
t_list = []
i=0.01
th_spar = 0.001
th_web = 0.002
i_list = []
while i <= 51.73/2:
    area = enclosed_area(i)
    polar = polar_moment_of_inertia(i,area,th_spar, th_web)
    torsionstiffnes = polar*G
    t_list.append(torsionstiffnes)
    i_list.append(i)
    i = i+0.1


#print(t_list)

# plt.plot(i_list,t_list, label = " polar inertia")
# plt.ylabel("m^4")
# plt.xlabel("span")
# plt.legend()
# plt.show()

# print(enclosed_area(10))
# print(z1(10),z4(10),x3(10),x1(10))
# print(polar_moment_of_inertia(10,enclosed_area(10),2.7))
