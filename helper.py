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
    def __init__(self, f:float, hl:list, Q:list, Z:list, A:list, binding:list ) -> None:
        """
        INPUT
        ---------
        f : hit frequency [Hz]
        hl : list of half life's [s]
        Q : list of Q values [keV]
        Z : list of Z values of daughter [1]
        A : list of A values of daughter [1]
        binding : list of values for the binding energy of daughter [keV]
        """
        hl = np.asarray(hl)
        Q = np.asarray(Q)
        Z = np.asarray(Z)
        A = np.asarray(A)
        binding = np.asarray(binding)

        self.f = f
        self.hl = hl
        self.Q = Q/1000 * EV_TO_J
        self.z = 2
        self.Z = Z
        self.A = A
        self.binding = binding/1000 * EV_TO_J
    
    def bc(self, Q, Z):
        return Q / (K* self.z * Z)
    
    def m(self, A, binding):
        mx = A*U_TO_KG + binding/C**2
        return MASS_ALPHA*mx / (MASS_ALPHA+mx)
    
    def G(self, Q, Z, binding, A, a):
        x = a / self.bc(A, binding)

        c1 = np.sqrt(2 * self.m(A, binding) / HBAR**2 / Q)
        c2 = self.z * Z * K
        c3 = np.arccos(np.sqrt(x)) - np.sqrt(x*(1-x))

        return c1*c2*c3
    
    def P(self, Q, Z, binding, A, a):
        return np.exp(-2*self.G(Q, Z, binding, A, a))
    
    def fit_func(self, x, a):
        Q, Z, binding, A = x
        return LOG2 / (self.f * self.P(Q, Z, binding, A, a))
    
    def fit(self):
        x = (self.Q, self.Z, self.binding, self.A)
        popt, pcov = optimize.curve_fit(self.fit_func, x, self.hl, bounds=([0],[1]))
        print(popt)
        print(pcov)
