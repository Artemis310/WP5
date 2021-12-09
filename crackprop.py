import Buckling as Bk
import numpy as np

class DamageTolerance:
    def __init__(self):
        self.K_IC = 29e6
        self.K_I = self.K_IC
        self.a = 5e-3

    def stress_allowed(self):
        return self.K_I/np.sqrt(np.pi*self
        )