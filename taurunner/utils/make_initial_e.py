import os
import numpy as np

from typing import Union

from .sample_powerlaw import sample_powerlaw

def make_initial_e(
    nevents: int,
    energy: Union[float, str],
    e_min=None,
    e_max=None
):
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
    Returns
    -------
    output (array-like) : Energies in eV
    '''
    # Make injected energies
    if isinstance(energy, str):
        if not os.path.isfile(energy):
            raise RuntimeError(f"Spline file {energy} does not exist")
        # sample initial energies and incoming angles from GZK parameterization
        if energy[-4:]=='.npy':
            cdf = np.load(energy, allow_pickle=True, encoding='latin1').item()
        elif energy[-4:]=='.pkl':
            import pickle as pkl
            with open(energy, 'rb') as pkl_f:
                cdf = pkl.load(pkl_f)
        cdf_indices = np.random.uniform(low=0., high=1.,size=nevents)
        eini = cdf(cdf_indices)
    else:
        e = float(energy)
        if e<=0:
            eini = sample_powerlaw(e_min, e_max, e + 1, size=nevents)
        else:
            eini = np.full(nevents, e)
    return eini
