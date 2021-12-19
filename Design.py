'''File used to configure and iterate the design of the wingbox'''

import Buckling as Bk
import TwistDistribution_requires_editing as Tw
import crackprop as Cp
import MMOI as Mm
import numpy as np
import matplotlib.pyplot as plt
import colorama as Cl
import moment_inerta as Mi

# Wingbox Constants
len_wbox = 51.73/2

# Set component properties

str_height = 0.25 # Stringer height [m]
str_width = 0.25  # Stringer width [m]
str_thick = 0.015 # Stringer thickness [m]
str_area = (str_height + str_width) * str_thick # Stringer area [m^2]


skin_thick_top = 0.01 # Top skin thickness [m]
skin_thick_bot = 0.01 # Bottom skin thickness [m]

flg_th = (skin_thick_top + skin_thick_bot) / 2 # Flange thickness [m]
spr_th = 0.012 # Spar thickness [m]

centroid_x, centroid_y = str_width / 2, str_height / 2 # Cross section centroids [m]

# Defining functions to be used in Design Configuring stage

def list_input(inputted_list):
    converted_to_list = inputted_list.split(" ")
    return list(map(float, converted_to_list))
            

# ------------------------------------------------- Design Configuring ---------------------------------------------------------------

#rib_loc = np.array(list_input(input("Input list of rib locations (make sure to include 0 and 25.865):"))) # Set locations of the ribs. Space separate list

rib_loc = [0, 3, 6, 9.05, 13, 15, 20, 25.865]

n_bays = len(rib_loc) - 1 # Number of bays

#top_str_bay_count = np.array(list_input(input("Input list of stringer count per bay (top):"))) # Number of stringers in each bay, make sure len(top_Str_bay) = n_bays
#bot_str_bay_count = np.array(list_input(input("Input list of stringer count per bay (bottom):"))) # Number of stringers in each bay, make sure len(bot_str_bay) = n_bays

top_str_bay_count = [12, 10, 10, 8, 6, 4, 3]
bot_str_bay_count = [10, 8, 8, 6, 4, 3, 3]

top_width_bay_str = np.zeros((1000, n_bays)) # Distances between top stringers in each bay, note each column is one bay
bot_width_bay_str = np.zeros((1000, n_bays)) # Same as line above but for bottom stringers

width_of_bay = np.zeros((1000, n_bays)) # Width of each bay, note each column is a bay

# This loop generates the distances between the stringers in each bay as well as the width of each bay

for i in range(n_bays):
    # POI = Bk.corner_points(top_str_bay_count[i], bot_str_bay_count[i], str_width, str_area, centroid_x, centroid_y, spr_th, flg_th,
    #                                     str_height, str_thick, rib_loc[i], rib_loc[i+1])
    span_values = np.linspace(rib_loc[i], rib_loc[i+1], 1000)
    
    bay_width = np.abs(Mi.x3(span_values))
    width_of_bay[:,i] = bay_width

    top_width_bay_str[:,i] = (bay_width - top_str_bay_count[i]*str_width)/(top_str_bay_count[i] - 1)
    bot_width_bay_str[:,i] = (bay_width - bot_str_bay_count[i]*str_width)/(bot_str_bay_count[i] - 1)

# Determination of kc ratios for top and bottom plates

max_plate_ratio_top = [] # Maxiumum a/b ratio for the top plate of each bay
max_plate_ratio_bot = [] # Maxiumum a/b ratio for the bottom plate of each bay

for i in range(n_bays):
    max_plate_ratio_top.append((rib_loc[i+1] - rib_loc[i])/min(top_width_bay_str[:,i]))
    max_plate_ratio_bot.append((rib_loc[i+1] - rib_loc[i])/min(bot_width_bay_str[:,i]))

print(f"Max a/b ratio of top plate per bay: {(max_plate_ratio_top)}")
print(f"Max a/b ratio of bot plate per bay: {(max_plate_ratio_bot)}")

# kc_top_bay = np.array(list_input(input("Input list of kc values (top):"))) # Input kc values for each bay, make sure len(kc_top_bay) = n_bays. Space separate the numbers
# kc_bot_bay = np.array(list_input(input("Input list of kc values (bot):"))) # Same as line above but for bottom bay

kc_top_bay = np.array([4, 4, 4, 4, 4, 4, 4])
kc_bot_bay = np.array([4, 4, 4, 4, 4, 4, 4])

# Determination of ks ratios for website

max_plate_ratio_web = [] # Maximum a/b ratio of the spar of each bay

for i in range(n_bays):
    span_values = np.linspace(rib_loc[i], rib_loc[i+1], 1000)
    
    max_plate_ratio_web.append(max((rib_loc[i+1] - rib_loc[i])/Mi.z1(span_values)))
    
print(f"Max a/b ratio of web per bay: {(max_plate_ratio_web)}")

# ks_per_bay = np.array(list_input(input("Input list of ks values:"))) # Input ks values for each bay, make sure len(kc_top_bay) = n_bays

ks_per_bay = np.array([5, 5, 5, 5, 5, 5, 5])

# k_str = int(input("Set the stringer buckling coefficient:")) # Stringer buckling coefficient

k_str = 1

# -------------------------------------------- Checking design for failure -------------------------------------------------------


# Crack propagation check
crck_check = Cp.CrackProp(top_str_bay_count[0], bot_str_bay_count[0], str_width, str_area, centroid_x, centroid_y, spr_th, flg_th,
                          str_height, str_thick)

if crck_check.check_crackprop_fail():
    print(f"Crack Propagation design {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
else:
    print(f"Crack Propagation design {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")
    
# Cycle throug each bay and check all failure criteria. If any fails, abort and print("fails"). Also, for each bay, calculate
# and plot Mos. failure criteria: buckleweb, buckleskin, bucklecolumn, plate tension.

for i in range(n_bays):
    # Multiple use vairables for each bay
    span_min = rib_loc[i]
    span_max = rib_loc[i+1]
    
    corner_coords = Bk.corner_points_vec(top_str_bay_count[i], bot_str_bay_count[i],
                                            str_width, str_area, centroid_x, centroid_y,
                                            spr_th, flg_th, str_height, str_thick, span_min, span_max)
    
    
    normal_buckle_stress_along_bay = Bk.NormalStressCalcs("Combined", top_str_bay_count[i], bot_str_bay_count[i],
                                            str_width, str_area, centroid_x, centroid_y,
                                            spr_th, flg_th, str_height, str_thick, corner_coords[1], corner_coords[-2]).stress_along_span(
                                            span_min, span_max)

    
    # Buckle check: Skin
    bckl_skin_top = Bk.BuckleSkin(kc_top_bay[i], skin_thick_top, top_width_bay_str[:,i])
    #bckl_skin_bot = Bk.BuckleSkin(kc_bot_bay[i], skin_thick_bot, bot_width_bay_str[:,i])
    
    
    if (normal_buckle_stress_along_bay[:,1] <= bckl_skin_top.crit_buckle_skin()).all():
        print(f"BAY {i+1}, Top Skin thickness (Buckling) {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
    else:
        print(f"BAY {i+1}, Top Skin thickness (Buckling) {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")
        
    # if not np.any(stress_along_bay[:,1] - bckl_skin_bot.crit_buckle_skin() <= 0):
    #     print(f"BAY {i+1}, Bottom Skin thickness (Buckling) {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
    # else:
    #     print(f"BAY {i+1}, Bottom Skin thickness (Buckling) {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")

    # Buckle check: stringers
    bckl_str = Bk.BuckleColumn(k_str, Mm.C_stringer(str_width, str_height, str_width, str_thick)[0], span_max - span_min, str_area)
    
    if (normal_buckle_stress_along_bay[:,1] <= bckl_str.crit_buckle_stringer()).all():
        print(f"BAY {i+1}, Stringer design {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
    else:
        print(f"BAY {i+1}, Stringer design {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")
        
    # Buckle check: web
    bckl_web = Bk.BuckleWeb(spr_th, span_min, span_max)
    
    if bckl_web.total_shear(ks_per_bay[i])[1]:
        print(f"BAY {i+1}, Spar thickness {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
    else:
        print(f"BAY {i+1}, Spar thickness {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")
    
    # Plate tension check
    
    tens_anal = Bk.Tension_analysis(top_str_bay_count[i], bot_str_bay_count[i],
                                            str_width, str_area, centroid_x, centroid_y,
                                            spr_th, flg_th, str_height, str_thick, 1000, 310e6, span_min, span_max)
    
    if tens_anal.check_for_failure():
        print(f"BAY {i+1}, Bottom Skin thickness (Tension) {Cl.Fore.GREEN} SUFFICIENT {Cl.Style.RESET_ALL}: \u2714")
    else:
        print(f"BAY {i+1}, Bottom Skin thickness (Tension) {Cl.Fore.RED} INSUFFICIENT {Cl.Style.RESET_ALL}: \u274c")
        
    # Margin of Safety plot, not including crack propagation    
    
#     mos = Bk.MarginOfSafety(top_str_bay_count[i], bot_str_bay_count[i],
#                                             str_width, str_area, centroid_x, centroid_y,
#                                             spr_th, flg_th, str_height, str_thick, span_min, span_max, k_str, span_max - span_min,
#                                             Mm.C_stringer(str_width, str_height, str_width, str_thick)[0], kc_top_bay[i],skin_thick_top,
#                                             width_of_bay[:, i])
    
#     span_values_plot = np.linspace(span_min, span_max, 1000)
#     plt.plot(span_values_plot, mos.plot_mos())
    
# plt.show()




class MarginOfSafety:
    def __init__(self, n_str_top, n_str_bot, str_width, str_height, str_thick, skin_thick_top, skin_thick_bot, spar_thick, rib_locations, top_str_count,
                 bot_str_count, flange_thick):
        self.n_str_top = n_str_top
        self.n_str_bot = n_str_bot 
        self.str_width = str_width
        self.str_height = str_height
        self.str_thick = str_thick
        self.skin_thick_top = skin_thick_top
        self.skin_thick_bot = skin_thick_bot
        self.spar_thick = spar_thick
        self.rib_locations = rib_locations
        self.top_str_count = top_str_count
        self.bot_str_count = bot_str_count
        self.str_area = (self.str_width, self.str_height) * self.str_thick
        self.centroid_x, self.centroid_y = self.str_width/2, self.str_height/2
        self.flange_thick = flange_thick
        
        
    def find_mos(self):
        for i in range(n_bays):
            norm_stress = Bk.NormalStressCalcs("Combined", self.n_str_top[i], self.n_str_bot[i], self.str_width, self.str_area,
                                               self.centroid_x, self.centroid_y, self.spar_thick, self.flange_thick, self.str_height,
                                               self.str_thick, Bk.corner_points[0], Bk.corner_points[-2]).stress_along_span(rib_loc[i], rib_loc[i+1])

            shear_stress = Bk.BuckleWeb(self.spar_thick, rib_loc[i], rib_loc[i+1]).total_shear()[0]
            
            skin_buckle = Bk.BuckleSkin(4, self.flange_thick, rib_loc[i+1] - rib_loc[i]).crit_buckle_skin()
            str_buckle = Bk.BuckleColumn(k_str, Mm.C_stringer(self.str_width, self.str_height, self.str_width, self.str_thick),
                                         rib_loc[i+1] - rib_loc[i], self.str_area).crit_buckle_stringer()
            web_buckle = Bk.BuckleWeb(self.spar_thick, rib_loc[i], rib_loc[i+1]).cri_buckle_web(ks_per_bay)
            
            
            