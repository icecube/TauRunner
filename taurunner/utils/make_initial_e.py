import os
import numpy as np
from taurunner.utils import units, sample_powerlaw, is_floatable

def make_initial_e(nevents, energy, e_min=None, e_max=None, rand=None):
    if rand is None:
        rand = np.random.RandomState()
    
    # Make injected energies
    if is_floatable(energy):
        e = float(energy)
        if e<=0:
            eini = sample_powerlaw(rand, e_min, e_max, e + 1, size=nevents)*units.GeV
        else:
            eini = np.full(nevents, e)*units.GeV
    else:
        # TODO figure out this part lol
        if not os.path.isfile(energy):
            raise RuntimeError(f"GZK CDF Spline file {energy} does not exist")
        # sample initial energies and incoming angles from GZK parameterization
        cdf_indices = rand.uniform(low=0., high=1.,size=nevents)
        cdf         = np.load(energy, allow_pickle=True, encoding='latin1').item()
        eini        = cdf(cdf_indices)*units.GeV
    return eini
