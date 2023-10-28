import numpy as np
from scipy import optimize


E_C = 1.602176634e-19
EPS = 8.8541878128e-12
HBAR = 1.054571817e-34
C = 299792458
MASS_ALPHA = 6.6446573357e-27
U_TO_KG = 1.66053906660e-27
EV_TO_J = 1.602176634e-19


K = E_C**2 / (4*np.pi*EPS)

LOG2 = np.log(2)

class Nuclei:
    def __init__(self, z, n, hl, hl_e, q, q_e, binding, binding_e):
        self.z = z
        self.n = n
        self.a = n+z
        self.hl = hl
        self.hl_e = hl_e
        self.q = q
        self.q_e = q_e
        self.binding = binding
        self.binding_e = binding_e
    
    def __lt__(self, other):
        return self.hl < other.hl
    
    def __eq__(self, other):
        return self.hl == other.hl
    
    def __gt__(self, other):
        return self.hl > other.hl

class Fit:
    def __init__(self, f:float, hl:float, Q:float, Z:float, A:float) -> None:
        self.f = f
        self.hl = hl
        self.Q = Q*1000 * EV_TO_J
        self.z = 2
        self.Z = Z
        self.A = A
        
        mx = A*U_TO_KG
        self.m = MASS_ALPHA*mx / (MASS_ALPHA+mx)

    def bc(self, Q, Z):
        return K * self.z * Z / Q
    
    def G(self, Q, Z, A, a):
        x = a / self.bc(Q, Z)

        c1 = np.sqrt(2 * self.m / HBAR**2 / Q)
        c2 = self.z * Z * K
        c3 = np.arccos(np.sqrt(x)) - np.sqrt(x*(1-x))

        return c1*c2*c3
    
    def P(self, Q, Z, A, a):
        return np.exp(-2*self.G(Q, Z, A, a))
    
    def fit_func(self, x, a):
        Q, Z = x
        return LOG2 / (self.f * self.P(Q, Z, self.A, a))
    
    def root_func(self, a):
        return self.fit_func((self.Q, self.Z), a) - self.hl
    
    def find_root(self):
        return optimize.root(self.root_func, 1.25e-15 * self.A**(1/3))