import numpy as np
from taurunner.modules import is_floatable
def make_initial_thetas(TR_specs, rand=None):

    if rand is None:
        rand=np.random.RandomState()

    if is_floatable(TR_specs['theta']):
        t = float(TR_specs['theta'])
        if t<0 or t>90:
            raise ValueError('Angles must be between 0 and 90')
        thetas = np.radians(np.ones(TR_specs['nevents'])*t)
    else:
        if TR_specs['theta']=='range':
            if TR_specs['th_min']<0 or TR_specs['th_max']>90:
                raise ValueError('Angles must be between 0 and 90')
            costh_min = np.cos(np.radians(TR_specs['th_max']))
            costh_max = np.cos(np.radians(TR_specs['th_min']))
            thetas = np.arccos(rand.uniform(costh_min, costh_max, TR_specs['nevents']))
        else:
            raise ValueError('theta sampling %s not suppoorted' % TR_specs['theta'])

    return thetas
