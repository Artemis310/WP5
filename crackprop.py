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

def allowed_stress(crack_length, KI):
    max_allowed_tension_stress = KI/((np.pi*crack_length)**(1/2))
    return max_allowed_tension_stress

max_allowed_stress = allowed_stress(0.005/2,29*10**6)

#print(f"max allowed stress = {max_allowed_stress/10e6:.2f} MPa")

span_stress_col = bk.Tension_analysis().stress_along_span()[-1]

#print(f"Actual tension stress at root = {span_stress_col[0,1]/10e6:.2f} MPa")

# if span_stress_col[0,1] > max_allowed_stress:
#     print("Crack prop req not met!")

# else:
#     print("You're all set!")

def check_crackprop_fail():
    if not max(span_stress_col) > max_allowed_stress:
        return True
    else:
        return False


#location = np.where(span_stress_col[:,1] > 1e6)

print(max(span_stress_col))
