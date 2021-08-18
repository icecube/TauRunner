import os
import numpy as np
from taurunner.utils import units, sample_powerlaw, is_floatable
def make_initial_e(TR_specs, rand=None):
    if rand is None:
        rand = np.random.RandomState()
    
    # Make injected energies
    if is_floatable(TR_specs['energy']):
        e = float(TR_specs['energy'])
        if e<=0:
            eini = sample_powerlaw(rand, TR_specs['e_min'], TR_specs['e_max'], e + 1, size=TR_specs['nevents'])*units.GeV
        else:
            eini = np.full(TR_specs['nevents'], e)*units.GeV
    else:
        if not os.path.isfile(TR_specs['energy']):
            raise RuntimeError(f"GZK CDF Spline file {TR_specs['energy']} does not exist")
        # sample initial energies and incoming angles from GZK parameterization
        cdf_indices = rand.uniform(low=0., high=1., size=TR_specs['nevents'])
        # TODO figure out this part lol
        cdf         = np.load(TR_specs['energy'], allow_pickle=True, encoding='latin1').item()
        eini        = cdf(cdf_indices)*units.GeV
    return eini
