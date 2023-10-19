"""

FIX UNITS

FIX MASS

"""


import numpy as np
from scipy import optimize

E_C = 1
EPS = 8.8541878128e-12
HBAR = 1
C = 1

K = E_C**2 / (4*np.pi*EPS)

LOG2 = np.log(2)

class Nuclei:
    def __init__(self, z, n, hl, hl_e, q, q_e):
        self.z = z
        self.n = n
        self.hl = hl
        self.hl_e = hl_e
        self.q = q
        self.q_e = q_e
    
    def __lt__(self, other):
        return self.hl < other.hl
    
    def __eq__(self, other):
        return self.hl == other.hl
    
    def __gt__(self, other):
        return self.hl > other.hl

class Fit:
    def __init__(self, f:float, hl:list, Q:list, Z:list, A:list, delta:list ) -> None:
        """
        f : hit frequency
        hl : list of half life's
        Q : list of Q values
        Z : list of Z values
        A : list of A values
        delta : list of values for the mass defect
        """
        hl = np.asarray(hl)
        Q = np.asarray(Q)
        Z = np.asarray(Z)
        A = np.asarray(A)
        delta = np.asarray(delta)

        self.f = f              # const
        self.hl = hl            # y
        self.Q = Q              # x1
        self.z = 2              # const
        self.Z2 = Z - 2         # x2
        self.A = A              # x3
        self.delta = delta      # x4
    
    def bc(self, Q, Z2):
        return Q / (K* self.z * Z2)
    
    def m(self, A, delta):
        return A + delta / C**2
    
    def G(self, Q, Z2, delta, A, a):
        x = a / self.bc(A, delta)

        c1 = np.sqrt(2 * self.m(A, delta) / HBAR**2 / Q)
        c2 = self.z * Z2 * K
        c3 = np.arccos(np.sqrt(x)) - np.sqrt(x*(1-x))

        return c1*c2*c3
    
    def P(self, Q, Z2, delta, A, a):
        return np.exp(-2*self.G(Q, Z2, delta, A, a))
    
    def fit_func(self, x, a):
        Q, Z2, delta, A = x
        return LOG2 / (self.f * self.P(Q, Z2, delta, A, a))
    

    def fit(self):
        x = (self.Q, self.Z2, self.delta, self.A)
        popt, pcov = optimize.curve_fit(self.fit_func, x, self.hl)
        print(popt)
        print(pcov)
