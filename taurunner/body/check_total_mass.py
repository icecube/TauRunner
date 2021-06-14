import numpy as np
from scipy.integrate import quad

def check_total_mass(Body):
    func = lambda x: 4*np.pi*Body.radius**3*Body.get_density(x)*x**2
    return quad(func, 0, 1, full_output=1)
