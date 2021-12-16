from scipy import interpolate
import numpy as np
import math
from MMOI import C_stringer
import matplotlib.pyplot as plt
import Shear_Calculations as scalc
from MMOI import C_stringer

# general constants
""" Here are the constants for the wingbox, these constants are related
to the geometry of the wingbox"""

c_spar1 = 0.2
c_spar2 = 0.61


# root chord formula
def c(y):
    """Function returning chord length as a function of the spanwise location in meters"""
    chord = 9.05 - (9.05 * (1 - 0.27)) / (51.73 / 2) * y
    return chord

c_vec =  np.vectorize(c)



# function coordinate top surface

""""Here we extract data from data files containing the coordinates of the required airfoil with chord 1. The data is 
then interpolated to create a continuous callable function"""

top_surf = open("top_surface.dat")

x1 = []
y1 = []

for line in top_surf.readlines():
    dummy = line.strip().split()
    y = float(dummy[1])
    x = float(dummy[0])
    x1.append(x)
    y1.append(y)

x1 = np.array(x1)
y1 = np.array(y1)

y_coord1 = interpolate.interp1d(x1, y1, kind="linear")

# function coordinate bottom surface

# function coordinate top surface

top_surf = open("bot_surface.dat")

x2 = []
y2 = []

for line in top_surf.readlines():
    dummy = line.strip().split()
    y = float(dummy[1])
    x = float(dummy[0])
    x2.append(x)
    y2.append(y)

x2 = np.array(x2)
y2 = np.array(y2)

y_coord2 = interpolate.interp1d(x1, y2, kind="linear")

# t/c constant
"""t/c ratios at the spar locations"""
t_c_spar1 = y_coord1(c_spar1) - y_coord2(c_spar1)
t_c_spar2 = y_coord1(c_spar2) - y_coord2(c_spar2)



def z1(y):
    dummy = t_c_spar1 * c(y)
    return dummy


def z4(y):
    dummy = t_c_spar2 * c(y)
    return dummy


def x1(y):
    dummy = (((c_spar2 - c_spar1) * c(y)) ** 2 + ((y_coord1(c_spar1) - y_coord1(c_spar2)) * c(y)) ** 2) ** 0.5
    return dummy


def x3(y):
    dummy = (c_spar2 - c_spar1) * c(y)
    return dummy


def spacing_str_top(y, n_str_top, width_str_top):
    dummy2 = (x1(y) - width_str_top * n_str_top) / (n_str_top - 1)
    return dummy2 #lmao Damien can't code


def spacing_str_bot(y, n_str_bot, width_str_bot):
    dummy2 = (x3(y) - width_str_bot * n_str_bot) / (n_str_bot - 1)
    return dummy2


# global scope constant for moment of inertia function
slope_top = math.atan((y_coord1(c_spar2) - y_coord1(c_spar1)) / (c_spar2 - c_spar1))


def moment_inertia_xx_func(y, n_str_top, n_str_bot, width_str, area_str, centroid_x, centroid_y, th_spar, th_flang, height_str,thick):
    # create array with verticale distance to centroid of stringer
    """th1 is the thickness for the spar and th2 for the flanges """
    # constants

    # initialise
    ds = spacing_str_top(y, n_str_top, width_str) + width_str
    d = y_coord1(c_spar1) * c(y) + centroid_x * math.sin(slope_top) - centroid_y * math.cos(slope_top)
    dd = ds * math.sin(slope_top)
    dist_arr = []

    while len(dist_arr) < n_str_top:
        dist_arr.append(d)
        d += dd
    dist_arr = np.array(dist_arr)

    """ Here under all the discrete/composite terms are computed for the neutral axis """
    centroid_steiner_terms_top = sum(dist_arr * area_str)
    centroid_term_x1 = x1(y) * th_flang * (0.5 * x1(y) * math.sin(slope_top) + y_coord1(c_spar1) * c(y))
    centroid_term_x3 = x3(y) * th_flang * y_coord2(c_spar1) * c(y)
    centroid_term_z1 = z1(y) * th_spar * (y_coord1(c_spar1) - y_coord2(c_spar1)) * c(y)
    centroid_term_z4 = z4(y) * th_spar * (y_coord1(c_spar2) - y_coord2(c_spar2)) * c(y)
    centroid_steiner_terms_bot = n_str_bot * area_str * (y_coord2(c_spar1) * c(y) + centroid_y)
    tot_area = (n_str_top + n_str_bot) * area_str + ((z1(y) + z4(y))) * th_spar + (x1(y) + x3(y)) * th_flang

    neutral_axis_height = (centroid_steiner_terms_top + centroid_steiner_terms_bot + centroid_term_z4 + centroid_term_z1 + centroid_term_x3 + centroid_term_x1) / tot_area

    """Now the neutral axis has been computed the moment of inertia terms will have to be computed"""
    moment_inert_string =  (n_str_top + n_str_bot) * C_stringer(width_str,height_str,width_str,thick)[0]
    inert_steiner_terms_top = sum((dist_arr - neutral_axis_height) ** 2 * area_str)
    inert_steiner_terms_bot = (y_coord2(c_spar1) * c(y) - neutral_axis_height) ** 2 * area_str * n_str_bot
    inert_term_x1 = (th_flang * x1(y) ** 3 * math.sin(slope_top) ** 2) / 12 + x1(y) * th_spar * (
            (0.5 * x1(y) * math.sin(slope_top) + y_coord1(c_spar1) * c(y)) - neutral_axis_height) ** 2
    inert_term_x3 = (th_flang * x3(y) ** 3) / 12 + x3(y) * th_spar * (y_coord2(c_spar1) * c(y) - neutral_axis_height) ** 2
    inert_term_z1 = (th_spar * z1(y) ** 3) / 12 + th_spar * (
            (y_coord1(c_spar1) * c(y) - y_coord2(c_spar1) * c(y)) - neutral_axis_height) ** 2
    inert_term_z4 = (th_spar * z4(y) ** 3) / 12 + th_spar * (
            (y_coord1(c_spar2) * c(y) - y_coord2(c_spar2) * c(y)) - neutral_axis_height) ** 2

    tot_inert = inert_steiner_terms_top + inert_steiner_terms_bot + inert_term_x3 + inert_term_z4 + inert_term_z1 + inert_term_x1 + moment_inert_string

    return neutral_axis_height, tot_inert


""""See example below, all relevant parameters are typed in for user friendly application"""

# inertia_xx = moment_inertia_xx_func(y= scalc.b_half , n_str_top=2, n_str_bot=2, width_str=0.1, area_str= 6 * 10 ** -4, centroid_x=0.025,
#                                     centroid_y=0.025, thickness=0.002)



def tot_area(y, n_str_top, n_str_bot, area_str, thickness):
    tot_area = (n_str_top + n_str_bot) * area_str + (z1(y) + z4(y) + x1(y) + x3(y)) * thickness
    return tot_area


def moment_inertia_yy_func(y, n_str_top, n_str_bot, width_str, area_str, centroid_x, centroid_y, th_spar, th_flang, height_str, thick):
    # create array with verticale distance to centroid of stringer
    """th1 is the thickness for the spar and th2 for the flanges """
    # constants


    """"The tip of the airfoil is used as datum"""
    # initialise top parameters
    ds_top = spacing_str_top(y, n_str_top, width_str) + width_str
    d = c_spar1 * c(y) +  centroid_x * math.cos(slope_top) + centroid_y * math.sin(slope_top)
    dd = ds_top * math.cos(slope_top)
    dist_arr_top = []

    # initialise bot parameters
    ds_bot = spacing_str_bot(y, n_str_bot, width_str) + width_str
    x = c_spar1 * c(y) + centroid_x
    dist_arr_bot = []

    while len(dist_arr_top) < n_str_top:
        dist_arr_top.append(d)
        d += dd
    dist_arr_top = np.array(dist_arr_top)

    while len(dist_arr_bot) < n_str_bot:
        dist_arr_bot.append(x)
        x += ds_bot
    dist_arr_bot = np.array(dist_arr_bot)

    """Again the neutral axis has to be computed  to enforce force equilibrium, we compute every composite term"""
    centroid_steiner_terms_top = sum(area_str * dist_arr_top)
    centroid_steiner_terms_bot = sum(area_str * dist_arr_bot)
    centroid_term_x1 = th_flang * x1(y) * (c_spar1 * c(y) + 0.5 * x1(y) * math.cos(slope_top))
    centroid_term_x3 = th_flang * x3(y) * (c(y) * c_spar1 + 0.5 * x3(y))
    centroid_term_z1 = th_spar * z1(y) * c_spar1 * c(y)
    centroid_term_z4 = th_spar * z4(y) * c_spar2 * c(y)
    tot_area = (n_str_top + n_str_bot) * area_str + ((z1(y) + z4(y))) * th_spar + (x1(y) + x3(y)) * th_flang

    neut_axis_x = (centroid_steiner_terms_top + centroid_steiner_terms_bot + centroid_term_z4 + centroid_term_z1 + centroid_term_x3 + centroid_term_x1) / tot_area

    moment_inert_string = (n_str_top + n_str_bot) * C_stringer(width_str, height_str, width_str, thick)[1]
    inert_steiner_terms_top = sum(area_str * (dist_arr_top - neut_axis_x) ** 2)
    inert_steiner_terms_bot =  sum(area_str * (dist_arr_bot - neut_axis_x) ** 2)
    inert_term_x1 = (th_flang * x1(y) ** 3 * math.cos(slope_top) ** 2) / 12 + x1(y) * th_spar * (neut_axis_x - (c_spar1 * c(y) + 0.5 * x1(y) * math.cos(slope_top))) ** 2
    inert_term_x3 = (x3(y) ** 3 * th_flang) / 12 + x3(y) * th_spar * ((c_spar1 * c(y) + 0.5 * x3(y)) - neut_axis_x) ** 2
    inert_term_z1 = (th_spar ** 3 * z1(y)) / 12 + z1(y) * th_spar * (neut_axis_x - c_spar1 * c(y)) ** 2
    inert_term_z4 = (th_spar ** 3 * z4(y)) / 12 + z4(y) * th_spar * (neut_axis_x - c_spar2 * c(y)) ** 2

    tot_inert = inert_steiner_terms_top + inert_steiner_terms_bot + inert_term_x3 + inert_term_z4 + inert_term_z1 + inert_term_x1 + moment_inert_string

    return neut_axis_x , tot_inert

#vectorisation
xx_vec_func = np.vectorize(moment_inertia_xx_func)
yy_vec_func = np.vectorize(moment_inertia_yy_func)
# span1 = np.linspace(0,25,50)
#
# plt.plot(span1, vec_func(span1, n_str_top=6, n_str_bot=4, width_str=0.1, area_str=0.5, centroid_x=0.03,
#                                    centroid_y=0.02, thickness=0.005)[1])
# plt.show()




