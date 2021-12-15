'''File used to configure and iterate the design of the wingbox'''

import Buckling as Bk
import TwistDistribution_requires_editing as Tw
import crackprop as Cp
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# Wingbox Constants
len_wbox = 51.73/2

# Set component properties

str_height = 1 # Stringer height [m]
str_width = 0.02  # Stringer width [m]
str_thick = 1 # Stringer thickness [m]
str_area = (str_height + str_width) * str_thick # Stringer area [m^2]

flg_th = 1 # Flange thickness [m]
spr_th = 1 # Spar thickness [m]

centroid_x, centroid_y = str_width / 2, str_height / 2 # Cross section centroids [m]


## Design configuring

rib_loc = np.array([0,8,len_wbox]) # Set locations of the ribs

n_bays = len(rib_loc) - 1 # Number of bays

top_str_bay_count = np.array([4,4]) # Number of stringers in each bay, make sure len(top_Str_bay) = n_bays
bot_str_bay_count = np.array([4,4]) # Number of stringers in each bay, make sure len(bot_str_bay) = n_bays


top_width_bay_str = np.zeros((1000, n_bays)) # Distances between top stringers in each bay, note each column is one bay
bot_width_bay_str = np.zeros((1000, n_bays)) # Same as line above but for bottom stringers


# This loop generates the distances between the stringers in each bay

for i in range(n_bays):
    POI = Bk.corner_points(top_str_bay_count[i], bot_str_bay_count[i], str_width, str_area, centroid_x, centroid_y, spr_th, flg_th,
                                        str_height, str_thick, rib_loc[i], rib_loc[i+1])

    bay_width = np.abs(POI[-1] - POI[-2])
    print(bay_width)

    top_width_bay_str[:,i] = (bay_width-top_str_bay_count[i]*str_width)/(top_str_bay_count[i] - 1)
    bot_width_bay_str[:,i] = (bay_width-bot_str_bay_count[i]*str_width)/(bot_str_bay_count[i] - 1)

# Determination of kc ratios for top and bottom plates

max_plate_ratio_top = []
max_plate_ratio_bot = []

for i in range(n_bays):
    max_plate_ratio_top.append((rib_loc[i+1] - rib_loc[i])/min(top_width_bay_str[:,i]))
    max_plate_ratio_bot.append((rib_loc[i+1] - rib_loc[i])/min(bot_width_bay_str[:,i]))

print(f"Max a/b ratio of top plate per bay: {max_plate_ratio_top}")
print(f"Max a/b ratio of bot plate per bay: {max_plate_ratio_bot}")

kc_top_bay = np.array([]) # Input kc values for each bay, make sure len(kc_top_bay) = n_bays
kc_bot_bay = np.array([]) # Same as line above but for bottom bay

# Determination of ks ratios for website

max_plate_ratio_web = []

'''for i in range(n_bays):
    print("bruh")'''


# Checking failure in each bay

# Cycle throug each bay and check all failure criteria. If any fails, abort and print("fails"). Also, for each bay, calculate
# and plot Mos. failure criteria: buckleweb, buckleskin, bucklecolumn, crackprop, plate tension.

for i in range(n_bays):
    # Multiple use vairables for each bay
    span_locations = np.linspace(rib_loc[i], rib_loc[i+1], 1000)
    span_min = rib_loc[i]
    span_max = rib_loc[i+1]
    
    corner_coords = Bk.corner_points_vec(top_str_bay_count[i], bot_str_bay_count[i],
                                            str_width, str_area, centroid_x, centroid_y,
                                            spr_th, flg_th, str_height, str_thick, span_min, span_max)
    
    
    stress_along_bay = Bk.NormalStressCalcs("Combined", top_str_bay_count[i], bot_str_bay_count[i],
                                            str_width, str_area, centroid_x, centroid_y,
                                            spr_th, flg_th, str_height, str_thick, corner_coords[0], corner_coords[-2]).stress_along_span(
                                            span_min, span_max)
    
    # Buckle check: Skin
    
    
    
    