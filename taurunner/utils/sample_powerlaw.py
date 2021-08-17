import numpy as np
def sample_powerlaw(rand, a, b, g, size=1):
    #Random spectrum function. g is gamma+1 (use -1 for E^-2)
    r = rand.uniform(low=0., high=1.0, size=size)
    if g == 0.:
        # E^-1 is uniform sampling in log space
        log_es = (np.log10(b) - np.log10(a)) * r + np.log10(a)
        return 10.**log_es
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)
