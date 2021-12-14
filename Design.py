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

n_bays = 2 # Set number of bays, equally spaced
len_bay = len_wbox / n_bays # Length of each bay [m]

rib_loc = np.array([0,8,len_wbox]) # Set locations of the ribs, make sure that len(rib_loc) = n_bays - 1

# This loop generates an array of all the rib locations
# for i in range(n_bays):
#     rib_loc = np.append(rib_loc, len_bay * (i+1))
    
print(rib_loc)
top_str_bay_count = np.array([4,4]) # Number of stringers in each bay, make sure len(top_Str_bay) = n_bays
bot_str_bay_count = np.array([4,4]) # Number of stringers in each bay, make sure len(bot_str_bay) = n_bays


top_width_bay_str = np.zeros((100, n_bays)) # Distances between top stringers in each bay, note each column is one bay
bot_width_bay_str = np.zeros((100, n_bays)) # Same as line above but for bottom stringers


# This loop generates the distances between the stringers in each bay

for i in range(n_bays):
    POI = Bk.corner_points(top_str_bay_count[i], bot_str_bay_count[i], str_width, str_area, centroid_x, centroid_y, spr_th, flg_th,
                                        str_height, str_thick, rib_loc[i], rib_loc[i+1])

    bay_width = np.abs(POI[-1] - POI[-2])
    print(bay_width)

    top_width_bay_str[:,i] = (bay_width-top_str_bay_count[i]*str_width)/(top_str_bay_count[i] - 1)
    bot_width_bay_str[:,i] = (bay_width-bot_str_bay_count[i]*str_width)/(bot_str_bay_count[i] - 1)

# Checking failure in each bay

