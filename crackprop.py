import Buckling as bk
import numpy as np
import seaborn as sns

# class DamageTolerance:
#     def __init__(self):
#         self.K_IC = 29e6
#         self.K_I = self.K_IC
#         self.a = 5e-3

#     def stress_allowed(self):
#         return self.K_I/np.sqrt(np.pi*self.a)



# #Bk.BuckleWeb.total_shear
# def half_crack_length(stress_applied,KI): #KI is crack thoughness -> material property KI=29*10**6
#     crack_length = ((KI/stress_applied)**2)/np.pi
#     total_crack_length = crack_length*2
#     return total_crack_length
class CrackProp:
    def __init__(self, n_str_top = 2, n_str_bot = 2, width_str = 0.02, 
                          area_str = 6e-4, centroid_x = 0.01, centroid_y = 0.01,
                          th_spar = 0.002, th_flang = 0.001, height_str = 0.02, thick = 0.002):
        self.n_str_top = n_str_top
        self.n_str_bot = n_str_bot
        self.width_str = width_str
        self.area_str = area_str
        self.centroid_x = centroid_x
        self.centroid_y = centroid_y
        self.th_spar = th_spar
        self.th_flang = th_flang
        self.height_str = height_str
        self.thick = thick
        self.span_stress_col = bk.Tension_analysis(self.n_str_top, self.n_str_bot, self.width_str, self.area_str, self.centroid_x, 
            self.centroid_y, self.th_spar, self.th_flang, self.height_str, self.thick).stress_along_span()[-1]
    
    def allowed_stress(crack_length, KI):
        max_allowed_tension_stress = KI/((np.pi*crack_length)**(1/2))
        return max_allowed_tension_stress

    max_allowed_stress = allowed_stress(0.005/2,29*10**6)
    
    def check_crackprop_fail(self):
        if not max(self.span_stress_col()) > self.max_allowed_stress:
            return True
        else:
            return False


#print(f"max allowed stress = {max_allowed_stress/10e6:.2f} MPa")



#print(f"Actual tension stress at root = {span_stress_col[0,1]/10e6:.2f} MPa")

# if span_stress_col[0,1] > max_allowed_stress:
#     print("Crack prop req not met!")

# else:
#     print("You're all set!")



#location = np.where(span_stress_col[:,1] > 1e6)

print(max(span_stress_col))
