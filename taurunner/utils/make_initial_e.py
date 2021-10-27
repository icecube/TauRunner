import os
import numpy as np
from taurunner.utils import units, sample_powerlaw, is_floatable

def make_initial_e(nevents, energy, e_min=None, e_max=None, rand=None):
    r'''
    Creates an array of energies
    Params
    ------
    nevents (int): Number of events
    energy (float or str): If positive float, then energy in GeV. If
        negative float, then spectral index. If str, then path
        to a flux cdf
    e_min (float or None): If not None, minimum energy in GeV 
        for a power-law
    e_max (float or None): If not None, maximum energy in GeV 
        for a power-law
    rand (np.random.RandomState): numpy random number generator object
    Returns
    -------
    output (array-like) : Energies in eV
    '''
    if rand is None:
        rand = np.random.RandomState()
    
    # Make injected energies
    if is_floatable(energy):
        e = float(energy)
        if e<=0:
            eini = sample_powerlaw(rand, e_min, e_max, e + 1, size=nevents)
        else:
            eini = np.full(nevents, e)
    else:
        if not os.path.isfile(energy):
            raise RuntimeError(f"Spline file {energy} does not exist")
        # sample initial energies and incoming angles from GZK parameterization
        if energy[-4:]=='.npy':
            cdf = np.load(energy, allow_pickle=True, encoding='latin1').item()
        elif energy[-4:]=='.pkl':
            import pickle as pkl
            with open(energy, 'rb') as pkl_f:
                cdf = pkl.load(pkl_f)
        cdf_indices = rand.uniform(low=0., high=1.,size=nevents)
        eini        = cdf(cdf_indices)
    return eini
