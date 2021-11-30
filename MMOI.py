import math
from math import pi


def L_stringer(height, lower_width, thickness):  # horizontal part of L is width and vertical part is height
    x_lower = lower_width / 2
    x_height = 0
    y_lower = 0
    y_height = height / 2

    x_bar = (lower_width * thickness * x_lower + height * thickness * x_height) / (
                lower_width * thickness + height * thickness)  # left bottom corner is datum feature
    y_bar = (lower_width * thickness * y_lower + height * thickness * y_height) / (
                lower_width * thickness + height * thickness)  # left bottom corner is datum feature

    tot_area = thickness * lower_width + thickness * height

    Ixx = (thickness * height ** 3) / 12 + height * thickness * (y_bar - y_height) ** 2 + lower_width * thickness * (
                y_bar - y_lower) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + lower_width * thickness * (
                x_bar - x_lower) ** 2 + height * thickness * (x_bar - x_height) ** 2

    return Ixx, Iyy, tot_area, x_bar, y_bar,


print("L-stringer: Ixx, Iyy, L_area, x_bar_L, y_bar_L = ", L_stringer(2, 1, 0.001))


def I_stringer(height, lower_width, upper_width, thickness):
    x_lower = 0
    x_height = 0
    x_upper = 0
    y_lower = 0
    y_height = height / 2
    y_upper = height

    x_bar = (lower_width * thickness * x_lower + height * thickness * x_height + upper_width * thickness * x_upper) / (
                lower_width * thickness + height * thickness + upper_width * thickness)  # bottom center is datum feature, a.k.a symmetric around y-axis
    y_bar = (height * thickness * y_height + lower_width * thickness * y_lower + upper_width * thickness * y_upper) / (
                height * thickness + lower_width * thickness + upper_width * thickness)  # bottom center is datum feature, a.k.a symmetric around y-axis

    tot_area = height * thickness + lower_width * thickness + upper_width * thickness

    Ixx = (thickness * height ** 3) / 12 + lower_width * thickness * (y_bar - y_lower) ** 2 + height * thickness * (
                y_height - y_bar) ** 2 + upper_width * thickness * (y_upper - y_bar) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + (thickness * upper_width ** 3) / 12
    return Ixx, Iyy, tot_area, x_bar, y_bar


print("I-stringer: Ixx, Iyy, I_area, x_bar_I, y_bar_I = ", I_stringer(2, 1, 1, 0.001))


def S_stringer(height, lower_width, width_upper, thickness):
    x_lower = -lower_width / 2
    x_height = 0
    x_upper = width_upper / 2
    y_lower = 0
    y_height = height / 2
    y_upper = height

    x_bar = (lower_width * thickness * x_lower + height * thickness * x_height + width_upper * thickness * x_upper) / (
            lower_width * thickness + height * thickness + width_upper * thickness)  # bottom center is datum feature
    y_bar = (lower_width * thickness * y_lower + height * thickness * y_height + width_upper * thickness * y_upper) / (
            lower_width * thickness + height * thickness + width_upper * thickness)  # bottom center is datum feature

    tot_area = lower_width * thickness + height * thickness + width_upper * thickness

    Ixx = (thickness * height ** 3) / 12 + lower_width * thickness * (y_bar - y_lower) ** 2 + height * thickness * (
                y_bar - y_height) ** 2 + width_upper * thickness * (y_bar - y_upper) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + (thickness * width_upper ** 3) / 12 + lower_width * thickness * (
                x_bar - x_lower) ** 2 + height * thickness * (x_bar - x_height) ** 2 + width_upper * thickness * (
                      x_bar - x_upper) ** 2
    return Ixx, Iyy, tot_area, x_bar, y_bar


print("S-stringer: Ixx, Iyy, S_area, x_bar_S, y_bar_S = ", S_stringer(2, 1, 1, 0.001))


def C_stringer(lower_width, height, upper_width, thickness):
    x_lower = lower_width / 2
    x_height = 0
    x_upper = upper_width / 2
    y_lower = 0
    y_height = height / 2
    y_upper = height

    x_bar = (lower_width * thickness * x_lower + height * thickness * x_height + upper_width * thickness * x_upper) / (
            lower_width * thickness + height * thickness + upper_width * thickness)  # left bottom corner is datum feature
    y_bar = (lower_width * thickness * y_lower + height * thickness * y_height + upper_width * thickness * y_upper) / (
            lower_width * thickness + height * thickness + upper_width * thickness)  # left bottom corner is datum feature

    tot_area = lower_width * thickness + height * thickness + upper_width * thickness

    Ixx = (thickness * height ** 3) / 12 + lower_width * thickness * (y_bar - y_lower) ** 2 + height * thickness * (
                y_bar - y_height) ** 2 + upper_width * thickness * (y_bar - y_upper) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + (thickness * upper_width ** 3) / 12 + lower_width * thickness * (
                x_bar - x_lower) ** 2 + height * thickness * (x_bar - x_height) ** 2 + upper_width * thickness * (
                      x_bar - x_upper) ** 2
    return Ixx, Iyy, tot_area, x_bar, y_bar


print("C-stringer: Ixx, Iyy, C_area, x_bar_C, y_bar_C = ", C_stringer(1, 1, 1, 0.001))


def T_stringer(height, lower_width, thickness):
    x_lower = 0
    x_height = 0
    y_lower = 0
    y_height = height / 2

    x_bar = (lower_width * thickness * x_lower + height * thickness * x_height) / (
                lower_width * thickness + height * thickness)  # center bottom is datum feature
    y_bar = (lower_width * thickness * y_lower + height * thickness * y_height) / (
                lower_width * thickness + height * thickness)  # center bottom is datum feature

    tot_area = thickness * lower_width + thickness * height

    Ixx = (thickness * height ** 3) / 12 + height * thickness * (y_bar - y_height) ** 2 + lower_width * thickness * (
                y_bar - y_lower) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + lower_width * thickness * (
                x_bar - x_lower) ** 2 + height * thickness * (x_bar - x_height) ** 2
    return Ixx, Iyy, tot_area, x_bar, y_bar


print("T-stringer: Ixx, Iyy, T_area, x_bar_T, y_bar_T = ", T_stringer(1, 1, 0.001))


def J_stringer(height, lower_width, radius_top, thickness):  # J stringer with round head
    x_lower = 0
    x_height = 0
    x_round = radius_top
    y_lower = 0
    y_height = height / 2
    y_round = height + (4 * radius_top) / (3 * math.pi)

    x_bar = (
                        lower_width * thickness * x_lower + height * thickness * x_height + math.pi * radius_top * thickness * x_round) / (
                        lower_width * thickness + height * thickness + math.pi * radius_top * thickness)  # center bottom is datum feature
    y_bar = (
                        lower_width * thickness * y_lower + height * thickness * y_height + math.pi * radius_top * thickness * y_round) / (
                        lower_width * thickness + height * thickness + math.pi * radius_top * thickness)  # center bottom is datum feature

    tot_area = lower_width * thickness + height * thickness + math.pi * radius_top * thickness

    Ixx = (thickness * height ** 3) / 12 + (1 / 8) * math.pi * radius_top ** 4 + lower_width * thickness * (
                y_bar - y_lower) ** 2 + height * thickness * (
                      y_bar - y_height) ** 2 + math.pi * radius_top * thickness * (y_bar - y_round) ** 2
    Iyy = (thickness * lower_width ** 3) / 12 + (1 / 8) * math.pi * radius_top ** 4 + lower_width * thickness * (
                x_bar - x_lower) ** 2 + height * thickness * (
                      x_bar - x_height) ** 2 + math.pi * radius_top * thickness * (x_bar - x_round) ** 2
    return Ixx, Iyy, tot_area, x_bar, y_bar


print("J-stringer: Ixx, Iyy, J_area, x_bar_J, y_bar_J = ", J_stringer(1, 1, 1, 0.001))
