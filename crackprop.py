#import Buckling as Bk
import numpy as np

class DamageTolerance:
    def __init__(self, K_IC, K_I, a = 5*10**(-3)):
        self.K_IC = K_IC
        self.K_I = K_I
        self.a = a

    def stress_allowed(self):
        return self.K_I/np.sqrt(np.pi*self.a)

#Bk.BuckleWeb.total_shear
def half_crack_length(stress_applied,KI): #KI is crack thoughness -> material property KI=29*10**6
    crack_length = ((KI/stress_applied)**2)/np.pi
    total_crack_length = a*2

def allowed_stress(crack_length, KI):
    max_allowed_tension_stress = KI/((np.pi*crack_length)**(1/2))
    return max_allowed_tension_stress





#def stress_crack_tip():
#test changes

