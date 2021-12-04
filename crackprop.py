import Buckling as Bk
import numpy as np

class DamageTolerance:
    def __init__(self, K_IC, K_I, a = 5*10**(-3)):
        self.K_IC = K_IC
        self.K_I = K_I
        self.a = a

    def stress_allowed(self):
        return self.K_I/np.sqrt(np.pi*self.a)
