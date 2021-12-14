'''File used to configure and iterate the design of the wingbox'''

import Buckling as Bk
import TwistDistribution_requires_editing as Tw
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# Wingbox Constants
len_wbox = 51.73/2

# Set component properties

str_height = 1 # Stringer height [m]
str_width = 1 # Stringer width [m]
str_thick = 1 # Stringer thickness [m]
str_area = (str_height + str_width) * str_thick # Stringer area [m^2]

flg_th = 1 # Flange thickness [m]
spr_th = 1 # Spar thickness [m]

centroid_x, centroid_y = 1,1 # Cross section centroids [m]


## Design configuring

n_bays = 2 # Set number of bays, equally spaced
len_bay = len_wbox / n_bays # Length of each bay [m]

rib_loc = np.array([0])

# This loop generates an array of all the rib locations
for i in range(n_bays):
    rib_loc = np.append(rib_loc, len_bay * (i+1))
    
print(rib_loc)
top_str_bay = np.array([2,2]) # Number of stringers in each bay, make sure len(top_Str_bay) = n_bays
bot_str_bay = np.array([2,2]) # Number of stringers in each bay, make sure len(bot_str_bay) = n_bays


top_width_bay_str = np.array([]) # Distances between top stringers in each bay, note each row is one bay
bot_width_bay_str = np.array([]) # Same as line above but for bottom stringers



# This loop generates the distances between the stringers in each bay
for i in range(n_bays):
    POI = Bk.corner_points(top_str_bay[i], bot_str_bay[i], str_width, str_area, centroid_x, centroid_y, spr_th, flg_th,
                                        str_height, str_thick, rib_loc[i], rib_loc[i+1])

    bay_width = np.abs(POI[-1] - POI[-2])
    
    if len(top_width_bay_str) == 0:
        top_width_bay_str = bay_width/(top_str_bay[i])
        bot_width_bay_str = bay_width/(bot_str_bay[i])
    else:
        top_width_bay_str = np.concatenate((top_width_bay_str, bay_width/(top_str_bay[i])), axis = 0)
        bot_width_bay_str = np.concatenate((bot_width_bay_str, bay_width/(bot_str_bay[i])), axis = 0)
    
print(top_width_bay_str)
