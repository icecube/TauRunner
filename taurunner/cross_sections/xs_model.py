import numpy as np
import pickle as pkl
from importlib.resources import path as imppath

import taurunner as tr
from taurunner.utils import units

NUTYPES      = ['nu', 'nubar']
NUCLEONS     = ['n', 'p']
INTERACTIONS = ['CC', 'NC']
TOT_DIF      = ['dsde', 'sigma']


# Helper functions
def tot_xs(E, spl):
    return np.exp(spl(np.log(E)))

def diff_xs(E_in, zz, spl):
    #E_min = 1e9 # Lowest knot on spline in eV
    #zz    = (E_out-E_min)/(E_in-E_min)
    res   = np.exp(spl(np.log(E_in), zz)[0])/E_in
    return res

class XSModel(object):

    def __init__(self, model:str, path:str=''):
        
        self.model = model
        if not path:
            with imppath('taurunner.resources.cross_section_tables', '__init__.py') as p:
                path = str(p).split('__init__.py')[0]


        # Iterate through all posibilities
        for nucleon in NUCLEONS:
            for td in TOT_DIF:
                for nutype in NUTYPES:
                    for interaction in INTERACTIONS:
                        desc_str = f'{nutype}_{nucleon}_{td}_{interaction}'
                        # TODO add check to throw more readable error
                        with open(f'{path}/{model}_{desc_str}.pkl', 'rb') as pkl_f:
                            setattr(self, f'_{desc_str}', pkl.load(pkl_f))

    def total_cross_section(self, E, nutype, interaction, proton_fraction=0.5):
        neutron_fraction = 1.0-proton_fraction
        val = proton_fraction *tot_xs(E, getattr(self, f'_{nutype}_p_sigma_{interaction}')) + \
              neutron_fraction*tot_xs(E, getattr(self, f'_{nutype}_n_sigma_{interaction}'))
        return val

    def differential_cross_section(self, Ein, zz, nutype, interaction, proton_fraction=0.5):
        zz = np.atleast_1d(zz)
        neutron_fraction = 1.0-proton_fraction
        val = proton_fraction *diff_xs(Ein, zz, getattr(self, f'_{nutype}_p_dsde_{interaction}')) + \
              neutron_fraction*diff_xs(Ein, zz, getattr(self, f'_{nutype}_n_dsde_{interaction}'))
        #val[Eout>Ein] = 0.
        if val.shape == (1,):
            val = val[0]
        return val
    
